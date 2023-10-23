import pandas as pd 
import os,re
from tqdm import tqdm

DATA_PATH = "/Users/shaderein/Library/CloudStorage/OneDrive-Personal/HKU/BDD"

exp = pd.DataFrame()

vehicle_path = os.path.join(DATA_PATH, "vehicle_experiment", "3Veh_ExpTask_Formal")

# human experiment

human_path = os.path.join(DATA_PATH, "human_experiment", "3Hum_ExpTask_Formal")
vehicle_path = os.path.join(DATA_PATH, "vehicle_experiment", "3Veh_ExpTask_Formal")

for root_path in [human_path, vehicle_path]:
    print("\n"+root_path.replace(DATA_PATH,''))
    test_dirs = [x for x in os.listdir(root_path) if os.path.isdir(os.path.join(root_path,x))]
    for test_dir in tqdm(test_dirs):
        test_num = int(re.findall(r"\d",test_dir)[0])
        test_dir = os.path.join(root_path,test_dir)

        # TODO: confirm with subject foler naming with subejct_no

        path = os.path.join(test_dir,"results")
        subject_dirs = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x)) and x!="Test"]

        assert(len(subject_dirs)==60)
        for subject_dir in tqdm(subject_dirs):
            subject_num = int(re.findall(r"\d+",subject_dir)[0])
            subject_dir = os.path.join(path,subject_dir)

            result = pd.read_csv(os.path.join(subject_dir,"RESULTS_FILE.txt"),delim_whitespace=True)

            result["test_num"] = [test_num for i in range(result.shape[0])]
            
            # TODO: check category used with Alice
            if "category" not in result:
                if "human" in root_path:
                    result["category"] = ["human" for i in range(result.shape[0])]
                elif "vehicle" in root_path:
                    result["category"] = ["car" for i in range(result.shape[0])]
                else:
                    assert(False)

            exp = pd.concat([exp,result])
            # exp.head()

    # exp.head()

    exp.head()
    exp.to_csv("exp_all_bdd.csv")