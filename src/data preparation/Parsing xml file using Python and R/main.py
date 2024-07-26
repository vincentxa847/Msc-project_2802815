import os
import xml.etree.ElementTree as ET


## Fetch the xml file of each experiment
# only for ebipipline
def file_path() :
    xmls = []
    directory = "./ebiPipeline/ebiPipeline/"
    srrs = os.listdir(directory)
    print(len(srrs)) # 23
    for srr in srrs :
        xml =  directory + srr + "/analysisConfig.xml"
        if os.path.exists(xml) :
            # print("File exists")
            xmls.append(xml)
        else :
            print(xml + "File not exists") # pfal3D7_invitro_versus_mosquito_produced_sporozoite_RNASeq don't have xml file, but fine it not using srr

    return xmls
xmls = file_path() # 22


## Fetch the SRR in xml file
def srr():
    srrs = []
    for xml in xmls:
        tree = ET.parse(xml)
        root = tree.getroot()  # Initialize the root
        for srr in root.iter("value"):
            srrs.append(srr.text)
    return srrs

SRR = srr()

## For the replicate, add the replicate number to make it unique
condition_counts = {} # condition is key and its count is value
final = []

for check in SRR:
    condition, unique_id = check.split("|") # "condition" (meaningful), not unique. "to replace" srr and err, unique
    if condition in condition_counts:
        condition_counts[condition] += 1
    else:
        condition_counts[condition] = 1

    # Create a unique condition string with replicate number if needed
    unique_condition = condition if condition_counts[condition] == 1 else f"{condition}__{condition_counts[condition]}"
    final.append(f"{unique_condition}|{unique_id}")
    

# Open a file in write mode
with open('SRRname.txt', 'w') as file:
    # Convert list elements to strings and join with newline characters
    lines = '\n'.join(map(str, final))
    # Write the lines to the file
    file.writelines(lines)
