import torch
import pandas as pd
from torch.utils.data import DataLoader, Dataset
from models import BCL_Network  # 假设这是你的自定义模型类

# 自定义数据集类
class MyDataSet(Dataset):
    def __init__(self, input):
        self.input = input

    def __len__(self):
        return len(self.input)

    def __getitem__(self, idx):
        x = torch.tensor(self.input.iloc[idx].values, dtype=torch.float32)
        x = x.unsqueeze(0)  # 增加一个通道维度，确保输入是 [channels, width]
        return x

# 预测并保存结果的函数
def predict_and_save(model_path, test_DataLoader, device, output_csv_path):
    # 加载模型
    model = BCL_Network().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()

    all_predictions = []
    with torch.no_grad():
        for x in test_DataLoader:
            x = x.to(device)
            output = model(x)
            predictions = output.cpu().numpy()
            all_predictions.extend(predictions)

    # 将预测结果保存为CSV
    results_df = pd.DataFrame({'Predicted Score': all_predictions})
    results_df.to_csv(output_csv_path, index=False)
    print(f"Predictions saved to {output_csv_path}")

# 运行预测的函数
def run_prediction():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    test_Batch_Size = 128  # 定义批处理大小

    # 数据文件路径
    test_data_path = "D:/Programming/python/PycharmProjects/ProteinDNABinding/data/rawdata/predict/test.txt"  # 替换为实际数据文件路径
    output_csv_path = "predictions.csv"  # 预测结果保存路径

    # 读取新的测试数据
    test_data = pd.read_csv(test_data_path, sep='\t')

    # 假设所有列都是特征
    X_test = test_data

    # 转换数据类型为数值类型
    X_test = X_test.apply(pd.to_numeric, errors='coerce')

    # 检查是否有NaN并处理NaN
    if X_test.isnull().values.any():
        print("Test data contains NaN values, filling with 0")
        X_test = X_test.fillna(0)

    # 打印输入数据的形状
    print(f"Shape of test data: {X_test.shape}")

    # 创建测试集数据加载器
    test_DataLoader = DataLoader(MyDataSet(X_test), batch_size=test_Batch_Size, shuffle=False)

    # 打印DataLoader的大小
    print(f"Test DataLoader size: {len(test_DataLoader.dataset)}")

    # 进行预测并保存结果
    best_model_path = "D:/Programming/python/PycharmProjects/ProteinDNABinding/model/NRF22/validate_params_5_4.pkl"
    predict_and_save(best_model_path, test_DataLoader, device, output_csv_path)

# 调用运行预测的函数
run_prediction()
