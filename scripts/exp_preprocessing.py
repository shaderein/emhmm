import pandas as pd 
import os,re
from tqdm import tqdm

DATA_PATH = "/mnt/h/OneDrive - The University Of Hong Kong/bdd/original/experiments_finished/"

exp = pd.DataFrame()

data_path = {
    "hum": {
        "data" : os.path.join(DATA_PATH, "human_experiment", "3Hum_ExpTask_Formal"),
        "block_start": [0,49,98,146]
        },
    "veh": {
        "data" : os.path.join(DATA_PATH, "vehicle_experiment", "3Veh_ExpTask_Formal"),
        "block_start": [0,49,98,148]
        }
}

for cat, data in data_path.items():
    root_path = data['data']
    block_start = data['block_start']

    print("\n"+root_path.replace(DATA_PATH,''))
    block_dirs = [x for x in os.listdir(root_path) if os.path.isdir(os.path.join(root_path,x))]
    for block_dir in tqdm(block_dirs):
        block_num = int(re.findall(r"\d",block_dir)[0])
        block_dir = os.path.join(root_path,block_dir)

        # TODO: confirm with subject foler naming with subejct_no

        path = os.path.join(block_dir,"results")
        subject_dirs = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x)) and x!="Test"]

        # assert(len(subject_dirs)==60)
        for subject_dir in tqdm(subject_dirs):
            subject_num = int(re.findall(r"\d+",subject_dir)[0])
            subject_dir = os.path.join(path,subject_dir)

            result = pd.read_csv(os.path.join(subject_dir,"RESULTS_FILE.txt"),sep='\t')

            result["block_num"] = [block_num for i in range(result.shape[0])]
            result['path'] = [subject_dir.replace(root_path,'') for i in range(result.shape[0])]

            # # Map trial to actual image id
            # #   the trial orders are appended as suffix in `results/saved_images` (0-based)
            # img_labels = os.listdir(os.path.join(subject_dir,"saved_images"))
            # img_labels_ord = [int(re.findall(r"_\d+", img)[0].replace("_","")) if re.findall(r"_\d+", img) else 0 for img in img_labels]
            # img_labels = [re.findall(r"^.+\.jpg", img)[0].replace("_hum",'') for img in img_labels]

            # img_labels_idx_sorted = sorted(range(len(img_labels_ord)), key=lambda k: img_labels_ord[k])
            # img_labels_sorted = [img_labels[i] for i in img_labels_idx_sorted]

            # result['image'] = img_labels_sorted

            # Edge case: incomplete files results in
            #   /vehicle_experiment/3Veh_ExpTask_Formal/Veh_ExpTask_Test4_bustruck_deploy/results/035ET4
            if not [a for a in os.listdir(os.path.join(subject_dir)) if re.findall(r"actual_TRIAL.+",a)]:
                print(f"\nIncomplete result:\t{subject_dir}\n")
                continue

            actual_trial_path = [a for a in os.listdir(os.path.join(subject_dir)) if re.findall(r"actual_TRIAL.+",a)][0]
            images_mapping = pd.read_csv(os.path.join(subject_dir,actual_trial_path),sep="\s+",header=None).iloc[:,2].tolist()
            result['image'] = [img.replace("_hum","") for img in images_mapping]

            result['trialID'] = [t + block_start[block_num-1] for t in result['Trial_Index_']]
            
            if "category" not in result:
                if cat == 'hum':
                    result["category"] = ["human" for i in range(result.shape[0])]
                elif cat == 'veh':
                    result["category"] = ["car" for i in range(result.shape[0])]
                else:
                    assert(False)

            exp = pd.concat([exp,result])
    exp.head()
    exp.to_csv(f"exp_{cat}_bdd.csv")