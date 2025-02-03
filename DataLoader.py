from datasets import load_dataset, Dataset
from typing import Dict, List

class MathProblemsDataset:
    def __init__(self, dataset_name : str, partition_name: str):
        self.data_partition = partition_name
        self.__ds = load_dataset(dataset_name, split=partition_name).remove_columns("messages")

    @property
    def data(self) -> Dataset:
        return self.__ds
    
    @data.setter #property-data.setter decorator
    def data(self, new_dataset: Dataset):
        self.__ds = new_dataset

    def get_size(self) -> int:
        return len(self.__ds)

    def get_problem(self, idx : int) -> str:
        return self.__ds['problem'][idx]

    def get_solution(self, idx : int) -> str:
        return self.__ds['solution'][idx]

    def get_problems(self) -> List[str]:
        return self.__ds['problem']

    def get_solutions(self) -> List[str]:
        return self.__ds['solution']
