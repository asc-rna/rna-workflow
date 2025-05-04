import pyarrow as pa
import pyarrow.ipc as ipc

def merge_arrow_files(file_list):
    # 创建一个空的 Table 列表
    tables = []

    # 读取每个文件并将其转换为 Table
    for file_name in file_list:
        with open(file_name, 'rb') as f:
            reader = ipc.RecordBatchReader(f)
            table = reader.read_all()
            tables.append(table)

    # 使用 concatenate_tables 函数合并所有 Table
    merged_table = pa.concat_tables(tables)

    return merged_table

def write_arrow_file(table, output_file):
    # 将合并后的 Table 写入一个新的 .arrow 文件
    with open(output_file, 'wb') as f:
        writer = ipc.RecordBatchFileWriter(f, table.schema)
        writer.write_table(table)
        writer.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-files", nargs=4, required=True, help="4 input .arrow files")
    parser.add_argument("-o", "--output-file", required=True, help="Output .arrow file")
    args = parser.parse_args()

    # 合并 .arrow 文件
    merged_table = merge_arrow_files(args.input_files)

    # 写入合并后的 .arrow 文件
    write_arrow_file(merged_table, args.output_file)