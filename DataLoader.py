from datasets import load_dataset
from typing import Dict, List


class MathDataset:
    def __init__(self, dataset_name : str, partition_name: str):
        self.data_partition = partition_name
        self.ds = load_dataset(dataset_name)[partition_name]
    
    def get_data(self) -> dict:
        return self.ds
    
    def get_problem(self, idx : int) -> str:
        return self.ds['problem'][idx]
    
    def get_solution(self, idx : int) -> str:
        return self.ds['solution'][idx]
    
    def get_problems(self) -> List[str]:
        return self.ds['problem']
    
    def get_solutions(self) -> List[str]:
        return self.ds['solution']
        
# Get data

ds = MathDataset(dataset_name="AI-MO/NuminaMath-TIR", partition_name='train')

dataset = ds.get_data()
problems = dataset['problem']

print(problems[:4])

# print("PROBLEM:")
# print(ds.get_problem(0))

# print("SOLUTION:")
# print(ds.get_solution(0))

