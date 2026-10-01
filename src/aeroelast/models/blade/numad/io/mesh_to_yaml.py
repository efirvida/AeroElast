import numpy as np
import yaml
from yaml import CLoader as Loader


def mesh_to_yaml(meshData, file_name):
    """
    TODO docstring
    """
    mDataOut = {}
    nodes = []
    for nd in meshData["nodes"]:
        ndstr = str(list(nd))
        nodes.append(ndstr)
    elements = []
    for el in meshData["elements"]:
        elstr = str(list(el))
        elements.append(elstr)
    esList = []
    for es in meshData["sets"]["element"]:
        newSet = {}
        newSet["name"] = es["name"]
        labels = []
        for el in es["labels"]:
            labels.append(int(el))
        newSet["labels"] = labels
        esList.append(newSet)
    nsList = []
    try:
        for ns in meshData["sets"]["node"]:
            newSet = {}
            newSet["name"] = ns["name"]
            labels = []
            for nd in ns["labels"]:
                labels.append(int(nd))
            newSet["labels"] = labels
            nsList.append(newSet)
    except Exception:
        pass
    sections = []
    for sec in meshData["sections"]:
        newSec = {}
        newSec["type"] = sec["type"]
        newSec["elementSet"] = sec["elementSet"]
        if sec["type"] == "shell":
            newLayup = []
            for lay in sec["layup"]:
                laystr = str(lay)
                newLayup.append(laystr)
            newSec["layup"] = newLayup
            newSec["xDir"] = str(list(sec["xDir"]))
            newSec["xyDir"] = str(list(sec["xyDir"]))
        else:
            newSec["material"] = sec["material"]
        sections.append(newSec)
    elOri = []
    for ori in meshData["elementOrientations"]:
        elOri.append(str(list(ori)))

    mDataOut["nodes"] = nodes
    mDataOut["elements"] = elements
    mDataOut["sets"] = {}
    mDataOut["sets"]["element"] = esList
    mDataOut["sections"] = sections
    mDataOut["elementOrientations"] = elOri

    try:
        mDataOut["materials"] = meshData["materials"]
    except Exception:
        pass

    fileStr = yaml.dump(mDataOut, sort_keys=False)

    fileStr = fileStr.replace("'", "")
    fileStr = fileStr.replace('"', "")

    outStream = open(file_name, "w")
    outStream.write(fileStr)
    outStream.close()


def yaml_to_mesh(fileName):
    inFile = open(fileName, "r")
    meshData = yaml.load(inFile, Loader=Loader)
    inFile.close()
    meshData["nodes"] = np.array(meshData["nodes"])
    meshData["elements"] = np.array(meshData["elements"])
    return meshData
