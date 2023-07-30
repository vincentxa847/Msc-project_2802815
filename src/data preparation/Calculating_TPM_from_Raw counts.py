#### VeuPathPipeline Raw count data into TPM (running TPMtool.py in python) 
import os
import subprocess

## Fetch the file name
def file_path(experiment) :
  directory = "./VeuPathPipeline/VeuPathPipeline/" + experiment
  directory = directory + "/samplesOutput/"
  sample = os.listdir(directory)
  updated_list = [directory + s for s in sample]
  return updated_list


# list of experiment
bartfai_time_series = file_path("Bartfai_time_series")

duffy = file_path("Duffy")

otto = file_path("Otto")

Su_seven_stages = file_path("Su_seven_stages")

Su_strand_specific = file_path("Su_strand_specific")

VEuPathPipline = [bartfai_time_series,duffy,otto,Su_seven_stages,Su_strand_specific]
VEuPathPipline = sum(VEuPathPipline, []) # unlist

TPM_output_names = []
for file_name in VEuPathPipline :
  names = file_name.split("/")
  TPM_output_name = "TPM." + names[7] + "." + names[9]
  TPM_output_names.append(TPM_output_name)


## Normalisation, using the script to transform raw count to TPM
for file_name in VEuPathPipline:
  # Deal with stranded data
  if "firststrand" in file_name:
    first_strand = file_name
    second_strand = file_name.replace("firststrand","secondstrand")

    names = file_name.split("/")
    TPM_output_name = "TPM." + names[7] + "." + names[9]
    first_strand_name = TPM_output_name
    second_strand_name = TPM_output_name.replace("firststrand","secondstrand")
    
    subprocess.run(['python', 'TPMtool.py','--genome','PlasmoDB-63_Pfalciparum3D7.gff','--input',first_strand,'--output',first_strand_name,'--stranded','--antisense_input',second_strand,'--antisense_output',second_strand_name], text=True)

  elif "secondstrand" in file_name:
    pass

  # Deal with unstranded data
  else:
    unstrand = file_name
    names = file_name.split("/")
    TPM_output_name = "TPM." + names[7] + "." + names[9]
    unstrand_name = TPM_output_name
        
    subprocess.run(['python', 'TPMtool.py','--genome','PlasmoDB-63_Pfalciparum3D7.gff','--input',unstrand,'--output',unstrand_name], text=True)


#### ebiPipeline Raw count data into TPM (running TPMtool.py in python) 

import os
import subprocess

## Normalisation, using the script to transform raw count to TPM
def file_path() :
    directory = "./Project_Data/ebiPipeline/ebiPipeline/"

    # using dataset and srr number to iterate
    experiments = os.listdir(directory)
    directory_dataset = [directory + experiment for experiment in experiments]
    directory_dataset_srr = []
    file_list = []

    for dataset in directory_dataset:
        srrs = os.listdir(dataset)
        # not using `os.path.join` here because it uses `\\` instead `/` to seperate
        directory_dataset_srr.append([dataset + "/" + srr for srr in srrs if os.path.isdir(os.path.join(dataset, srr))])

    for listofsample in directory_dataset_srr:
        for sample in listofsample:
            files = os.listdir(sample)
            # not using `os.path.join` here because it uses `\\` instead `/` to seperate
            file_list.append([sample + "/" + file for file in files])

    return file_list

ebiPipline = file_path()
ebiPipline = sum(ebiPipline, []) # unlist

## Normalisation, using the script to transform raw count to TPM
for file_name in ebiPipline:

    # Deal with stranded data
    if "firststrand" in file_name:
        first_strand = file_name
        second_strand = file_name.replace("firststrand", "secondstrand")

        names = file_name.split("/")
        TPM_output_name = names[7] + "." + names[9]
        first_strand_name = TPM_output_name
        first_strand_name = first_strand_name.replace("counts", "TPM")
        second_strand_name = first_strand_name.replace("firststrand", "secondstrand")

        subprocess.run(
            ['python', 'TPMtool.py', '--genome', 'PlasmoDB-63_Pfalciparum3D7.gff', '--input', first_strand, '--output',
             first_strand_name, '--stranded', '--antisense_input', second_strand, '--antisense_output',
             second_strand_name], text=True)

    elif "secondstrand" in file_name:
        pass

    # Deal with unstranded data
    else:
        unstrand = file_name
        names = file_name.split("/")
        TPM_output_name = names[7] + "." + names[9]
        TPM_output_name = TPM_output_name.replace("counts", "TPM")
        unstrand_name = TPM_output_name

        subprocess.run(
            ['python', 'TPMtool.py', '--genome', 'PlasmoDB-63_Pfalciparum3D7.gff', '--input', unstrand, '--output',
             unstrand_name], text=True)
