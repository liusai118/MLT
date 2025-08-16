import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, f1_score

# 假设 BCL_Network 已经定义在 models 模块中
from models import BCL_Network
import myDataSet as ms

# 准备数据加载器的函数
def getFullDataLoader(X, y, batch_size):
    full_DataSet = ms.MyDataSet(input=X.reset_index(drop=True), label=y.reset_index(drop=True))
    full_DataLoader = DataLoader(dataset=full_DataSet, batch_size=batch_size, shuffle=False)
    return full_DataLoader

# 预测函数
def predict(model_path, test_loader):
    # 加载模型
    model = BCL_Network().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()  # 切换到评估模式

    all_outputs = []

    with torch.no_grad():
        for inputs, _ in test_loader:
            inputs = inputs.to(device)
            outputs = model(inputs)
            # 将回归值转换为二进制标签
            binary_outputs = (outputs > 0.5).int()
            all_outputs.append(binary_outputs.cpu().numpy())

    # 确保预测结果是1维数组
    return np.concatenate(all_outputs).flatten()

# 保存预测结果为CSV
def save_predictions(sequences, y_pred, file_path):
    predictions_df = pd.DataFrame({
        "Sequence": sequences,
        "Predicted Label": y_pred
    })
    predictions_df.to_csv(file_path, index=False)
    print(f"Predictions saved to {file_path}")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# 加载数据
all_data = pd.read_csv('D:\\Programming\\python\\PycharmProjects\\ProteinDNABinding\\data\\rawdata\\NRF22\\test3', sep='\t')

# 打印数据的前几行和列名进行调试
print(all_data.head())
print(all_data.columns)

# 提取序列数据和标签
sequences = all_data.iloc[:, 0]
X = sequences
y = all_data.iloc[:, 1]

# 假设 batch_size 是之前定义好的、·
batch_size = 128
test_loader = getFullDataLoader(X, y, batch_size)

# 假设 model_path 是保存模型的路径
model_path = 'D:/Programming/python/PycharmProjects/ProteinDNABinding/NRF22/fold_3/validate_params_epoch_2.pkl'

# 进行预测
y_pred = predict(model_path, test_loader)

# 保存预测结果为CSV
output_csv_path = 'D:\\Programming\\python\\PycharmProjects\\ProteinDNABinding\\results\\predictions3.csv'
save_predictions(sequences, y_pred, output_csv_path)
