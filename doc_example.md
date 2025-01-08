<!-- 仅供写README的参考 -->
# SparTA

## Install
准备cuda12.1，创建一个新的python=3.10虚拟环境
```
source /home/zhaijidong/scy/spack/share/spack/setup-env.sh
spack load python@3.10
spack load cuda@12.1

python -m venv /path/to/venv/store/location
source </path/to/venv/store/location> + "/bin/activate"
```
启元集群上的cuda12.1在`/data/apps/tools/cuda/cuda-12.1/`目录下

注意spack最新版不会自动加载一些必要的环境变量，手动加载或查阅文档解决。

在python环境中用pip安装对应cu12.1的torch==2.5.0。

clone Sparta (https://github.com/zheng-ningxin/SparTA.git) 仓库，进入pit_sm80分支，使用pip编译安装
```
pip install -e . --verbose
```

注意，安装过程中可能同时需要网络和GPU，在启元集群上，可以salloc一个节点，并通过两层ssh-remote-forwarding让计算节点联网。

```
salloc --nodes=1 --ntasks=8 --time=23:00:00 --partition=arch --cpus-per-task=8 --gres=gpu:8

ssh <g3010> -R 10817:localhost:10817  //我本地的代理端口首先转发到了登陆节点的10817端口，再转发到计算节点的10817端口

port=10817
export https_proxy=http://127.0.0.1:$port 
export http_proxy=http://127.0.0.1:$port
export all_proxy=socks5://127.0.0.1:$port
```


测试
```
python ./test_sparse_moe_fp16.py
```
注意参照Note修改参数和真实专家路由张量的存储路径。

## Note
本仓库基于pit_sm80的修改：

1. 注释掉了setup.py中与sparse_moe无关的模块，减少修改cuda代码后重新编译安装的时间

2. 对应1.，注释了./sparta/opset/\_\_init__.py 中的部分（被删除的模块的）import代码。

3. ./csrc/moe_sparse_forward_kernel.cu，修改了FP16的forward_function函数的kernel启动代码；在此处修改inhidden和outhidden参数，与下一点中的测试脚本对齐（修改后需要重新编译）。

4. 从./test/test_sparse_moe_fp16.py中修改出了与我们对齐的测试脚本，在./test_sparse_moe_fp16.py位置。使用torch计算图，主要可调变量见下。
```py
    B = 1 #BatchSize
    S = 128 #不含Begin Token的Token lenth
    N_exp = 64 //专家数
    in_hidden =  2048 #MLP in size ，需与cuda代码对齐
    out_hidden = 176 #MLP middle size，需与cuda代码对齐
    topk = 6 #选取排名前k的专家进行路由转发
    exp_dis_path = "dump/inputlen128/router_logits_layer5_1.pt" #真实专家路由分数，inputLlen需要与上面的S(sequence len)一致
```



5. 参考`DynamicSparseMoE`定义了`SparseMoEAdvisor`类，与主流的LLM推理层对齐，包含pit手写的MLP矩阵乘法和silu激活，mul，shared exp等组件，结构见下面的代码。
```py
    def forward(self, tokens, expids):

        sparse_moe.convert_index(expids, self.sparse_index, self.expert_count)
        with torch.no_grad():
            x1 = sparse_moe.forward(tokens, self.weight, expids, self.sparse_index, self.expert_count,self.GLOBAL_M)
            x2 = sparse_moe.forward(tokens, self.weight2, expids, self.sparse_index, self.expert_count,self.GLOBAL_M)
            
            x1 = siluLayer(x1)
            x1 = torch.mul(x1,x2)
            
            moe_res = sparse_moe.forward(x1,self.reverse_weight,expids,self.sparse_index,self.expert_count,self.GLOBAL_M)

            shared_mid = self.shared1(tokens)
            shared_mid = siluLayer(shared_mid)
            shared_res = self.shared2(shared_mid)

            res = moe_res + shared_res
            # res = moe_res 
        return res
```