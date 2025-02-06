from datasets import load_dataset, Dataset
from typing import Dict, List

class MathProblemsDataset:
    def __init__(self, dataset_name: str, partition_name: str):
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
    
    def exemine_clusters(self, cluster_col: str, show_solution: bool = False) -> None:
        print(type(self.__ds[cluster_col][:5]))
        for cluster in set(self.__ds[cluster_col]):
            print(40*"#")
            print(f"======= Cluster {cluster} ==============")
            print(40*"#")
            for problem in self.__ds.filter(lambda example: example[cluster_col]==cluster).select([0,1,2,3,4]):
                print(3*'===============', 'PROBLEM:')
                print(problem['problem'])

                if show_solution:
                    print(3*'===============', 'SOLUTION:')
                    print(problem['solution'])
