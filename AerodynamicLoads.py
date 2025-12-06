import numpy as np

def read_CL_from_file(filename):
    with open(filename, "r") as f:
        for line in f:
            if line.strip().startswith("CL"):
                # Split on '=' and convert the right side to float
                return float(line.split('=')[1])
    raise ValueError("CL value not found in file.")

CL0 = read_CL_from_file("a0.txt")
CL10 = read_CL_from_file("a10.txt")

#A function that sorts the data from XFLR 5 to get the proper lists
def AL_sorter(data):
    y = []
    ch = []
    cl = []
    cd = []
    cm = []
    ai = []

    for i in range(len(data)):
        if data[i][0] >= 0:
            y.append(data[i][0].item())
            ch.append(data[i][1].item())
            cl.append(data[i][3].item())
            cd.append(data[i][5].item())
            cm.append(data[i][7].item())

    ylst = np.array(y)
    chlst = np.array(ch)
    cllst = np.array(cl)
    cdlst = np.array(cd)
    cmlst = np.array(cm)

    return ylst, chlst, cllst, cdlst, cmlst

#interpolates the data from XFLR 5 to any CL value
def AL_specific(f0, f10, CLL, CLL0 = CL0, CLL10 = CL10):

    slope = (f10 - f0) / (CLL10 - CLL0)
    cl_CL = f0 + slope * (CLL - CLL0)

    return cl_CL

#computes the angle of attack for a specific CL
def AL_computeaoa(CL):

    slope = (CL10 - CL0)/10
    aoa = (CL - CL0)/slope

    return aoa













