from datasets import load_dataset
from typing import Dict, List

class MathProblemsDataset:
    def __init__(self, dataset_name : str, partition_name: str):
        self.data_partition = partition_name
        self.ds = load_dataset(dataset_name, split=partition_name)

    def get_data(self) -> dict:
        return self.ds

    def get_size(self) -> int:
        return len(self.ds)

    def get_problem(self, idx : int) -> str:
        return self.ds['problem'][idx]

    def get_solution(self, idx : int) -> str:
        return self.ds['solution'][idx]

    def get_problems(self) -> List[str]:
        return self.ds['problem']

    def get_solutions(self) -> List[str]:
        return self.ds['solution']
