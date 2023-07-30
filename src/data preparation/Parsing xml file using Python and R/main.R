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
# replicate this process for three times
i = []
k = []
for check in SRR:
    two_string = check.split("|")
    condition = two_string[0]  # condition (meaningful), not unique
    to_replace = two_string[1]  # srr and err, unique
    if condition in i:
        add_mark_condition = condition + "__2"
        k.append(add_mark_condition + "|" +to_replace)
    else :
        k.append(condition + "|" + to_replace)

    i.append(condition)
j = []
z = []
for second_check in k :
    two_string = second_check.split("|")
    condition = two_string[0]  # condition (meaningful), not unique
    to_replace = two_string[1]  # srr and err, unique
    if condition in j:
        add_mark_condition = condition.replace("__2","__3")
        z.append(add_mark_condition + "|" + to_replace)
    else:
        z.append(condition + "|" + to_replace)
    j.append(condition)

n = []
final = []
for third_check in z :
    two_string = third_check.split("|")
    condition = two_string[0]  # condition (meaningful), not unique
    to_replace = two_string[1]  # srr and err, unique
    if condition in n:
        add_mark_condition = condition.replace("__3","__4")
        final.append(add_mark_condition + "|" + to_replace)
    else:
        final.append(condition + "|" + to_replace)
    n.append(condition)

# Open a file in write mode
with open('SRRname.txt', 'w') as file:
    # Convert list elements to strings and join with newline characters
    lines = '\n'.join(map(str, final))
    # Write the lines to the file
    file.writelines(lines)
