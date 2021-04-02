# **************************************************************
# Quantization and Visualization of Interfacial Chemistry in APT
# **************************************************************

# Given an Atom Probe Tomography (APT) dateset, this code quantifies 
# and visualizes the concentration of one element on the isosurface 
# of another element. Visualization in 3D requires package *mayavi*.

# Authors: 
#    Linqing Peng <linqingp@outlook.com> from Grinnell College
#    Jonathan D. Poplawsky <poplawskyjd@ornl.gov> at the Center for 
# Nanophase Materials Sciences, Oak Ridge National Laboratory (ORNL)

# Acknowledgements: 
#    This research was conducted at ORNL's Center for Nanophase 
# Materials Sciences, which is a DOE Office of Science User Facility. 
# L.P. participated in this work at ORNL under the Oak Ridge Science 
# Semester (ORSS) program. This research was supported in part by 
# ORNL's Laboratory Directed Research and Development program.

# Required input:

# 1. Range file .rrng.
# 2. Atom position file that contains the location of all of the 
#    atoms. It can either be a .csv file with x, y, z coordinates and 
#    mass-to-charge ratio of each atom, or a .wrl file that combines 
#    the x, y, z coordinates of each atom with the isosurface.
# 3. isosurface file .wrl. Can be combined with the atom positions.
# 4. The name of the element of interest.
# 5. The maximun distance of an atom of interest from the isosurface 
#    in nm.
# Please list the names of the input files or the inputs below.

# **************************************************************
# User Input
# **************************************************************

# Input name of range file
# Replace the string on the right side of the equal sign with the 
# full name of range file with extension
rangeFileName = "RR350 reconstruction.rrng"

# Input type of input
# Indicate by "YES" or "NO" if a .csv file is used for atom 
# positions (or they are included in the .wrl file for 
# isosurface)
useCSVInput = "NO"

# Input name of .csv file
# If you choose "YES" above, specify below the name of the .csv 
# file. If "NO", skip next line.
pointFileName = "R29_11916-v01.CSV"

# Input name of VRML file
# Replace the string on the right side of the equal sign with 
# the full name of VRML file with extensions
pointAndSurfaceFileName = "RR350_Cu_18perc_ISO_allPoints_1Vox.wrl"

# Input element of interest
# Replace the string on the right side of the equal sign with 
# the symbol of the element of interest, e.g. "Zr"
xElementName = "Mn"

# Input max distance of atoms from interface 
# Replace the float on the right side of the equal sign with 
# the maximun distance of an atom of interest from the 
# isosurface in nanometer
xDisNM = 1

# Input bin size to combine triangles when calcuating
# concentration
# Replace the integer on the right side of the equal sign with 
# how many times the bin size of xElement is larger than that 
# of the element used to create the isosurface
binSize = 4

# Input grid size used to store atoms based on position
# Indicate the voxel size in nanometer, normally don't need to 
# change
voxelSizeXNM = 1
voxelSizeYNM = 1
voxelSizeZNM = 1

# **************************************************************
# Preliminary processing of input
# **************************************************************

# Input the maximum distance requirement above and below the 
# isosurface for an atom to be counted. 
# In .wrl mode, the user input is in nanometer, so here we 
# convert the distance into the unit used in VRML file in this 
# case the conversion is 130.3 nm = "1"
# TO-DO: how to adapt to different project
if (useCSVInput.upper() == "YES"):
    useCSV = 1 # 1 represent "true"
else:
    useCSV = 0 # 0 represent "false"
    if (useCSVInput.upper() != "NO"):
        print("Cannot recognize useCSV input. Use default \"NO\""\
                "value instead")
        # We assume useCSVInput = "NO" here so that the program 
        # can continue even if users don't use correct input.

# Make sure xDisNM is non-negative
xDisNM = abs(xDisNM)

# Process range input file
# We assume the number of elements and that of ranges are both 
# nonnegative.
# We assume the lower and upper bounds of ranges are all positive
import re

# Open range input file
rangeInputObject = open(rangeFileName, "r")
rangeInput = rangeInputObject.read()

# Declare variables
elementNumber = 0 # the number of differnt types of elements 
                  # appeared in ions
elementName = [] # a list of strings that correspond to the names 
                 # of elements

# Input ion types
elementNumberMatch = re.search("\[Ions\](?:\n|\r\n)Number=(\d+)",\
        rangeInput, re.M) 
if elementNumberMatch is None:
    print("#Wrong number of elements.")
else:
    elementNumber = int(elementNumberMatch.group(1))

for i in range(0, elementNumber):
    elementNameMatch = re.search("Ion"+ str(i+1) + "=(\S+)", \
            rangeInput)
    if elementNameMatch is None:
        print("#Wrong name of elements.") 
    else:
        elementName.append(elementNameMatch.group(1))

# Declare variables
elementCombNumber = 0 # the number of different element	combinations 
                      # appeared in ions. Element combinations include
                      # molecules, so elementCombNumber = 
                      # elementNumber + number of molecules
rangeNumber = 0 # the number of different m/q ranges. rangeNumber = 
                # len(rangeMin) = len(rangeMax)
rangeMin = [] # a list of the lower m/q bounds of ranges
rangeMax = [] # a list of the upper m/q bounds of ranges
elementCombIndex = [] # a 2D list of the elements contained in each				
                      # combination. elementCombIndex[i][j] gives the 					
                      # index in elementName of j^th element in the 
                      # i^th combination
elementCombAbund = [] # a 2D list of the number of atoms of each element 
 				      # contained in one combination. elementCombAbund
                      # [i][j] gives the number of atoms of j^th element 
                      # contained in one i^th combination. For example, 
				      # one ion Al2O1 contains 2 Al.
rangeElementComb = [] # a list of the index of the element combination
				      # to which each range maps

# Input m/q range
rangeNumberMatch = re.search("\[Ranges\](?:\n|\r\n)Number=(\d+)", \
        rangeInput)
if rangeNumberMatch is None:
    print("#Wrong number of ranges.")
else:
    rangeNumber = int(rangeNumberMatch.group(1))

for i in range(0, rangeNumber):
    rangeMatch = re.search("Range" + str(i+1) + "=(\d+(\.\d*)?|\.\d+) "\
            "(\d+(\.\d*)?|\.\d+) Vol:(\d+(\.\d*)?|\.\d+) (.+) Color:", \
            rangeInput)
    if rangeMatch is None:
        print("#Wrong range #" + str(i+1) + ".") 
    else:
        rangeMin.append(float(rangeMatch.group(1)))
        rangeMax.append(float(rangeMatch.group(3)))
        combList = rangeMatch.group(7).split(' ')
        tempElementIndexList = []
        tempElementAbundList = []
        for j in range(len(combList)):
            componentMatch = re.search("(\S+):(\d+)", combList[j])
            found = 0
            for k in range(0, elementNumber):
                if componentMatch.group(1) == elementName[k]:
                    tempElementIndexList.append(k)
                    tempElementAbundList.append(int(\
                                componentMatch.group(2)))
                    found = 1
                    break
            if not(found):
                print("#Wrong element name for range " + str(i+1))
        found = 0
        for j in range(len(elementCombIndex)): 
            if (cmp(tempElementIndexList,elementCombIndex[j]) == 0 and \
                    cmp(tempElementAbundList,elementCombAbund[j]) == 0):
                rangeElementComb.append(j)
                found = 1
                break
        if not(found):
            rangeElementComb.append(len(elementCombIndex))
            elementCombIndex.append(tempElementIndexList)
            elementCombAbund.append(tempElementAbundList)
elementCombNumber = len(elementCombIndex)

# Close file to release space
rangeInputObject.close()

# Print debug info
#print("There are " + str(elementNumber) + " elements")
#print("Names of elements are " + str(elementName)) 
#print("There are " + str(elementCombNumber) + " element combinations")
#print("Combinations of elements are " + str(elementCombIndex))
#print("The number of atoms of each element contained in the \
#combinations are " + str(elementCombAbund))
#print("There are " + str(rangeNumber) + " ranges") 
#print("Lower bounds of ranges are " + str(rangeMin)) 
#print("Upper bounds of ranges are " + str(rangeMax)) 
#print("Elements ranges each correspond to are " + str(rangeElementComb)) 

# **************************************************************
# Process atoms and isosurface input file 
# **************************************************************

# WRL input atoms mode 
# All coordinates of atoms, vertices of isosurface, grid 
# information and xDis will be stored in normalized unit 
# so that all coordinates are in the range [-1, 1]


if not(useCSV):
    # open isosurface input file
    isosurfaceInputObject = open(pointAndSurfaceFileName, "r")
    isosurfaceInput = isosurfaceInputObject.read()

    # declare variables
    pointX = []
    pointY = []
    pointZ = [] # store the positions of atoms
    pointSetRegEx = re.compile("geometry PointSet")
    pointRegEx = re.compile("([-+]?(\d+(\.\d*)?|\.\d+)) "\
"([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    verticesRegEx = re.compile("([-+]?(\d+(\.\d*)?|\.\d+)) "\
"([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    coordIndexRegEx = re.compile("coordIndex")
    groupRegEx = re.compile("Group")

    # input points
    verticesStartPos = re.search("geometry IndexedFaceSet", \
            isosurfaceInput).end()
    for pointSetIterator in pointSetRegEx.finditer(isosurfaceInput, \
            0, verticesStartPos):
        pointX.append([])
        pointY.append([])
        pointZ.append([])
        pointSetEndPos = groupRegEx.search(isosurfaceInput, \
                pointSetIterator.end()).start()
        for pointIterator in pointRegEx.finditer(isosurfaceInput, \
                pointSetIterator.end(), pointSetEndPos): 
            pointX[-1].append(float(pointIterator.group(1)))
            pointY[-1].append(float(pointIterator.group(4)))
            pointZ[-1].append(float(pointIterator.group(7)))
    pointNumber = sum(len(pointX[i]) for i in range(elementCombNumber)) 

    # declare variables
    vertexX = []
    vertexY = [] 
    vertexZ = [] # store the positions of vertices of the 
  			# triangular grid

    # input vertices
    vertexTrRegEx = re.compile("translation ([-+]?(\d+(\.\d*)?"\
"|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    vertexTr = vertexTrRegEx.findall(isosurfaceInput, pointSetEndPos, \
            verticesStartPos)[-1]
    vertexTrX = float(vertexTr[0])
    vertexTrY = float(vertexTr[3])
    vertexTrZ = float(vertexTr[6])

    # find the scale, assuming the scale on all x, y, z directions are 
    # the same
    scaleStartPos = re.search("DEF Scales Switch", \
isosurfaceInput).end()
    gridSizeRegEx = re.compile("translation ([-+]?(\d+(\.\d*)?|"\
"\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    gridSize = gridSizeRegEx.search(isosurfaceInput, \
scaleStartPos)
    gridSizeX = 2*abs(float(gridSize.group(1)))
    gridSizeY = 2*abs(float(gridSize.group(4)))
    gridSizeZ = 2*abs(float(gridSize.group(7)))
    OriLenXRegEx = re.compile("string \[ \"X \(" + "([-+]?(\d+"\
"(\.\d*)?|\.\d+))")
    OriLenX = float(OriLenXRegEx.search(isosurfaceInput, \
scaleStartPos).group(1))
    OriLenYRegEx = re.compile("string \[ \"Y \(" + "([-+]?(\d+"\
"(\.\d*)?|\.\d+))")
    OriLenY = float(OriLenYRegEx.search(isosurfaceInput, \
scaleStartPos).group(1))
    OriLenZRegEx = re.compile("string \[ \"Z \(" + "([-+]?(\d+"\
"(\.\d*)?|\.\d+))")
    OriLenZ = float(OriLenZRegEx.search(isosurfaceInput, \
                scaleStartPos).group(1))
    if abs((gridSizeX) - 1.0) < 0.00000001:
        scale = abs(OriLenX / gridSizeX)
    elif abs((gridSizeY) - 1.0) < 0.00000001:
        scale = abs(OriLenY / gridSizeY)
    else:
        scale = abs(OriLenZ / gridSizeZ)
    gridMinX = gridSizeX / -2 
    gridMinY = gridSizeY / -2 
    gridMinZ = gridSizeZ / -2 

    #vertexScRegEx = re.compile("scale ([-+]?(\d+(\.\d*)?|\.+"\
        #"\d+)) ([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    #vertexSc = vertexScRegEx.findall(isosurfaceInput, \
#pointSetEndPos, verticesStartPos)[-1]
    #vertexScX = float(vertexSc[0])
    #vertexScY = float(vertexSc[3])
    #vertexScZ = float(vertexSc[6])

    triangleEndPos = re.search("normal Normal", isosurfaceInput).start()
    triangleStartPos = coordIndexRegEx.search(isosurfaceInput, \
            verticesStartPos, triangleEndPos).end()
    for verticesIterator in \
        verticesRegEx.finditer(isosurfaceInput, verticesStartPos, \
                triangleStartPos):
        vertexX.append(float(verticesIterator.group(1)) / \
                scale + vertexTrX) # replace "* vertexScX"
        vertexY.append(float(verticesIterator.group(4)) / \
                scale + vertexTrY) # replace "* vertexScY"
        vertexZ.append(float(verticesIterator.group(7)) / \
                scale + vertexTrZ) # replace "* vertexScZ"
    vertexNumber = len(vertexX)

    # Declare variables
    # store the connections between vertices to form triangles
    triangle = [] 

    # Input triangles
    triangleRegEx = re.compile("(\d+) (\d+) (\d+) -1")
    for triangleIterator in triangleRegEx.finditer(isosurfaceInput, \
            triangleStartPos, triangleEndPos):
        triangle.append([int(triangleIterator.group(1)), \
            int(triangleIterator.group(2)), int(triangleIterator.group(3))])
    triangleNumber = len(triangle)

    # Scale the xDis and voxel size to the proper unit
    xDis = abs(xDisNM / scale)
    voxelSizeX = abs(voxelSizeXNM / scale)
    voxelSizeY = abs(voxelSizeYNM / scale)
    voxelSizeZ = abs(voxelSizeZNM / scale)
    
    
    # Close file to release space
    isosurfaceInputObject.close()
    
    # Print debug info
    #print("There are " + str(pointNumber) + " points")
    #print("Arranged in catagories of elements: " + str(list(len(pointX[i])\
    #                 for i in range(elementCombNumber))))
    #print("X coordinates of atoms are " + str(pointX))
    #print("Y coordinates of atoms are " + str(pointY))
    #print("Z coordinates of atoms are " + str(pointZ))
    #print("There are " + str(vertexNumber) + " vertices")
    #print("First 10 X coordinates of vertices of the triangular "\
    #        "grid are " + str(vertexX[:10]))
    #print("Last 10 X coordinates of vertices of the triangular "\
    #        "grid are " + str(vertexX[-11:]))
    #print("First 10 Y coordinates of vertices of the triangular "\
    #        "grid are " + str(vertexY[:10]))
    #print("Last 10 Y coordinates of vertices of the triangular \
    #        grid are " + str(vertexY[-11:]))
    #print("First 10 Z coordinates of vertices of the triangular \
    #        grid are " + str(vertexZ[:10]))
    #print("Last 10 Z coordinates of vertices of the triangular \
    #        grid are " + str(vertexZ[-11:]))
    #print("There are " + str(triangleNumber) + " triangles")
    #print("First 10 connections of vertices of the triangular \
    #        grid are " + str(triangle[:10]))
    #print("Last 10 connection of vertices of the triangular \
    #        grid are " + str(triangle[-11:]))
    #print("Translations on x,y,z direction are " + str(vertexTrX) + "\
    #        "' ' + str(vertexTrY) + ' ' + str(vertexTrZ))
    #print("Scales on x,y,z direction are " + str(vertexScX) \
    #        + ' ' + str(vertexScY) + ' ' + str(vertexScZ))
    #print("One unit length in the wrl file is equal to " + \
    #        str(scale) + "nm in the sample")
    #print("Atoms within " + str(xDis) + " unit length above and below "\
    #    "the surface will be calculated")
    #print("The overall size of the rectangular grid is %.3f * %.3f * "\
    #        "%.3f" % (gridSizeX, gridSizeY, gridSizeZ))
    #print("The grid starts with %.3f * %.3f * %.3f" % \
    #        (gridMinX, gridMinY, gridMinZ))
    #print("The size of the rectangular voxel is %.6f * %.6f * %.6f" % \
    #        (voxelSizeX, voxelSizeY, voxelSizeZ))

else:
    # **************************************************************
    # Process atoms and isosurface input file 
    # **************************************************************

    # CSV input atoms mode
    # All coordinates of atoms, vertices of isosurface, grid information
    # and xDis will be stored in unit of nanometers

    # Iterate through all of the points in the csv file to input atom 
    # positions
    # Convert their mass-to-charge ratio to element type. 
    
    # Open atom position input file
    pointInputObject = open(pointFileName, "r")
    pointInput = pointInputObject.read()
    
    # Declare variables
    # Store the positions of atoms. Each sub-list embeded in the outer 
    # list includes all atoms of one element
    pointX = [[] for i in range(elementCombNumber)]
    pointY = [[] for i in range(elementCombNumber)]
    pointZ = [[] for i in range(elementCombNumber)] 
    #input points
    for pointIterator in re.finditer("([-+]?(\d+(\.\d*)?|\.\d+)),"\
            "([-+]?(\d+(\.\d*)?|\.\d+)),([-+]?(\d+(\.\d*)?|\.\d+)),"\
            "([-+]?(\d+(\.\d*)?|\.\d+))", pointInput):
        pointMtoE = float(pointIterator.group(10))
        for i in range(rangeNumber):
            if (pointMtoE > rangeMin[i]) and (pointMtoE < rangeMax[i]):
                pointX[rangeElementComb[i]].append(float(pointIterator.group(1)))
                pointY[rangeElementComb[i]].append(float(pointIterator.group(4)))
                pointZ[rangeElementComb[i]].append(float(pointIterator.group(7)))
                break
    pointNumber = sum(len(pointX[i]) for i in range(elementCombNumber))
    
    # Open isosurface input file
    isosurfaceInputObject = open(pointAndSurfaceFileName, "r")
    isosurfaceInput = isosurfaceInputObject.read()

    # Declare variables
    verticesRegEx = re.compile("([-+]?(\d+(\.\d*)?|\.\d+)) "\
            "([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    coordIndexRegEx = re.compile("coordIndex")
    groupRegEx = re.compile("Group")

    # Declare variables
    vertexX = []
    vertexY = []
    vertexZ = [] #store the positions of vertices of the triangular grid

    # Input vertices
    verticesStartPos = re.search("geometry IndexedFaceSet", \
            isosurfaceInput).end()
    triangleEndPos = re.search("normal Normal", isosurfaceInput).start()
    triangleStartPos = coordIndexRegEx.search(isosurfaceInput, \
            verticesStartPos, triangleEndPos).end()
    for verticesIterator in verticesRegEx.finditer(isosurfaceInput, \
            verticesStartPos, triangleStartPos):
        vertexX.append(float(verticesIterator.group(1))) 
        vertexY.append(float(verticesIterator.group(4)))
        vertexZ.append(float(verticesIterator.group(7)))
    vertexNumber = len(vertexX)

    # Declare variables
    triangle = [] #store the connections between vertices to form triangles

    # Input triangles
    triangleRegEx = re.compile("(\d+) (\d+) (\d+) -1")
    for triangleIterator in triangleRegEx.finditer(isosurfaceInput, \
            triangleStartPos, triangleEndPos):
        triangle.append([int(triangleIterator.group(1)), \
            int(triangleIterator.group(2)), int(triangleIterator.group(3))])
    triangleNumber = len(triangle)
    
    # Input the size of the grid     
    vertexTrRegEx = re.compile("translation ([-+]?(\d+(\.\d*)?|\.\d+)) "\
            "([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    vertexTr = vertexTrRegEx.findall(isosurfaceInput, 0, verticesStartPos)[-1]
    vertexTrX = float(vertexTr[0])
    vertexTrY = float(vertexTr[3])
    vertexTrZ = float(vertexTr[6])
    scaleStartPos = re.search("DEF Scales Switch", isosurfaceInput).end()
    gridSizeRegEx = re.compile("translation ([-+]?(\d+(\.\d*)?|\.\d+)) "\
            "([-+]?(\d+(\.\d*)?|\.\d+)) ([-+]?(\d+(\.\d*)?|\.\d+))")
    gridSize = gridSizeRegEx.search(isosurfaceInput, scaleStartPos)
    gridSizeX = 2*abs(float(gridSize.group(1)))
    gridSizeY = 2*abs(float(gridSize.group(4)))
    gridSizeZ = 2*abs(float(gridSize.group(7)))
    OriLenXRegEx = re.compile("string \[ \"X \(" + "([-+]?(\d+(\.\d*)?|\.\d+))")
    OriLenX = float(OriLenXRegEx.search(isosurfaceInput, scaleStartPos).group(1))
    OriLenYRegEx = re.compile("string \[ \"Y \(" + "([-+]?(\d+(\.\d*)?|\.\d+))")
    OriLenY = float(OriLenYRegEx.search(isosurfaceInput, scaleStartPos).group(1))
    OriLenZRegEx = re.compile("string \[ \"Z \(" + "([-+]?(\d+(\.\d*)?|\.\d+))")
    OriLenZ = float(OriLenZRegEx.search(isosurfaceInput, scaleStartPos).group(1))
    if abs((gridSizeX) - 1.0) < 0.00000001:
        scale = abs(OriLenX / gridSizeX)
    elif abs((gridSizeY) - 1.0) < 0.00000001:
        scale = abs(OriLenY / gridSizeY)
    else:
        scale = abs(OriLenZ / gridSizeZ)
    gridMinX = (gridSizeX / -2 - vertexTrX) * scale
    gridMinY = (gridSizeY / -2 - vertexTrY) * scale
    gridMinZ = (gridSizeZ / -2 - vertexTrZ) * scale
    gridSizeX = OriLenX
    gridSizeY = OriLenY
    gridSizeZ = OriLenZ
    
    

    # Scale the xDis and voxel size to the proper unit
    xDis = abs(xDisNM)
    voxelSizeX = abs(voxelSizeXNM)
    voxelSizeY = abs(voxelSizeYNM)
    voxelSizeZ = abs(voxelSizeZNM)

    # Close file to release space
    isosurfaceInputObject.close()
    pointInputObject.close()
    
    # Print debug info
    print("There are " + str(pointNumber) + " points")
    print("Arranged in catagories of elements: " + str(list(len(pointX[i]) \
                    for i in range(elementCombNumber))))
    #print("X coordinates of atoms are " + str(pointX))
    #print("Y coordinates of atoms are " + str(pointY))
    #print("Z coordinates of atoms are " + str(pointZ))
    print("There are " + str(vertexNumber) + " vertices")
    #print("First 10 X coordinates of vertices of the triangular grid are "\
    #        + str(vertexX[:10]))
    #print("Last 10 X coordinates of vertices of the triangular grid are " \
    #        + str(vertexX[-11:]))
    #print("First 10 Y coordinates of vertices of the triangular grid are " \
    #        + str(vertexY[:10]))
    #print("Last 10 Y coordinates of vertices of the triangular grid are " \
    #        + str(vertexY[-11:]))
    #print("First 10 Z coordinates of vertices of the triangular grid are " \
    #        + str(vertexZ[:10]))
    #print("Last 10 Z coordinates of vertices of the triangular grid are " \
    #        + str(vertexZ[-11:]))
    print("There are " + str(triangleNumber) + " triangles")
    #print("First 10 connections of vertices of the triangular grid are " \
    #        + str(triangle[:10]))
    #print("Last 10 connection of vertices of the triangular grid are " \
    #        + str(triangle[-11:]))
    #print("Translations on x,y,z direction are " + str(vertexTrX) + ' ' \
    #        + str(vertexTrY) + ' ' + str(vertexTrZ))
    #print("Scales on x,y,z direction are " + str(vertexScX) + ' ' \
    #        + str(vertexScY) + ' ' + str(vertexScZ))
    #print("The size of the rectangular grid is %.3f * %.3f * %.3f" % \
    #        (gridSizeX, gridSizeY, gridSizeZ))
    #print("The grid starts with %.3f * %.3f * %.3f" % \
    #        (gridMinX, gridMinY, gridMinZ))
    #print("The size of the rectangular voxel is %.6f * %.6f * %.6f" % \
    #        (voxelSizeX, voxelSizeY, voxelSizeZ))
    
# **************************************************************
# Map each atom into cubic grid. 
# **************************************************************

# Construct a hash table with separate chaining to resolve collision.
import math
import time

time0 = time.time()
# Define helper function for vector subtraction, dot product and cross 
# product in 3D
def subt(u,v): #u-v
    return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

def dotProd(u,v): #u dot v
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

def crossProd(u,v): #u cross v
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], \
                         u[0] * v[1] - u[1] * v[0]]

def det2D(u,v):
    return u[0] * v[1] - u[1] * v[0]

time1 = time.time()
print("Start calculation at " + str(time1 - time0))

# Decide how many voxels there are in the rectangular grid along x, 
# y and z directions
voxelNumX = int(math.ceil(gridSizeX / voxelSizeX))
voxelNumY = int(math.ceil(gridSizeY / voxelSizeY))
voxelNumZ = int(math.ceil(gridSizeZ / voxelSizeZ))
time2 = time.time()
print("Start to initialize 3D grid at " + str(time2 - time0))

# Declare variables
gridToPoint = [[[[[] for m in range(elementCombNumber)] for k in \
    range(voxelNumZ)] for j in range(voxelNumY)] for i in range(voxelNumX)]
time3 = time.time()
print("3D grid is initialized at " + str(time3 - time0))

for i in range(len(pointX)):
    for j in range(len(pointX[i])):
        gridToPoint[int(math.floor((pointX[i][j] - gridMinX) / voxelSizeX))]\
    [int(math.floor((pointY[i][j] - gridMinY) / voxelSizeY))]\
    [int(math.floor((pointZ[i][j] - gridMinZ) / voxelSizeZ))][i].append(j)
time4 = time.time()
print("3D grid is created at " + str(time4 - time0))
        
# Convert the symbol of the element that we will graph to its corresponding 
# index in the element combination list
xElementCombIndex = [] # a list of the index of the element combinations that 
                       # contain xElement named xElementName
xElementCombAbund = [] # a list of the number of atoms of xElement one 
                       # corresponding element oombination contains
for i in range(elementCombNumber):
    for j in range(len(elementCombIndex[i])):
        if xElementName == elementName[elementCombIndex[i][j]]:
            xElementCombIndex.append(i)
            xElementCombAbund.append(elementCombAbund[i][j])
            break

# Define variables
pointOfxElement = []
pointOfOtherElement = []
triangleValueOfxElement = []
triangleValueOfAllElement = []
triangleConcentrationOfxElement = []
isTriangle = []

# Define a "new" set of variables that exclude all the triangles 
# that have no area
triangleValueOfxElementNew = []
triangleValueOfAllElementNew = []
triangleConcentrationOfxElementNew = []
triangleNew = []


# Iterate through all triangles and for each triangle, look for points 
# in nearby voxels that fall into the triangular prism above&below the triangle
for i in range(triangleNumber):
    v = []
    for j in range(3):
        v.append([vertexX[triangle[i][j]], vertexY[triangle[i][j]], \
                vertexZ[triangle[i][j]]])
    v01 = subt(v[1], v[0])
    v02 = subt(v[2], v[0])
    n = crossProd(v01, v02) #normal vector of the triangle
    triangleValueOfxElement.append(0)
    triangleValueOfAllElement.append(0)
    # When two vertices have same coordinates or when three vertices are 
    # on the same line ("same" = indistinguishable due to float point 
    # precision), i.e. the three vertices don't define a plane, the normal 
    # vector is (0,0,0) and no point is above or below the triangle. Thus,
    # we will just skip this triangle.
    if dotProd(n,n) == 0:
        isTriangle.append(0)
        triangleConcentrationOfxElement.append(0)
        continue
    isTriangle.append(1)
    triangleNew.append(triangle[i])
    # Rectangular culling to narrow down the range of cubes we need 
    # to iterate through
    xMax = max(vertexX[triangle[i][j]] for j in range(3))
    yMax = max(vertexY[triangle[i][j]] for j in range(3))
    zMax = max(vertexZ[triangle[i][j]] for j in range(3))
    xMin = min(vertexX[triangle[i][j]] for j in range(3))
    yMin = min(vertexY[triangle[i][j]] for j in range(3))
    zMin = min(vertexZ[triangle[i][j]] for j in range(3))
    xVoxelMax = int(math.floor((xMax + xDis - gridMinX) / voxelSizeX + 1))
    yVoxelMax = int(math.floor((yMax + xDis - gridMinY) / voxelSizeY + 1))
    zVoxelMax = int(math.floor((zMax + xDis - gridMinZ) / voxelSizeZ + 1))
    xVoxelMin = int(math.floor((xMin - xDis - gridMinX) / voxelSizeX))
    yVoxelMin = int(math.floor((yMin - xDis - gridMinY) / voxelSizeY))
    zVoxelMin = int(math.floor((zMin - xDis - gridMinZ) / voxelSizeZ))
    for xVoxel in range(max(0, xVoxelMin), min(voxelNumX, xVoxelMax)):
        for yVoxel in range(max(0, yVoxelMin), min(voxelNumY, yVoxelMax)):
            for zVoxel in range(max(0, zVoxelMin), min(voxelNumZ, zVoxelMax)):
                for m in range(elementCombNumber):
                    # Decide whether points in the m^th range are xElement
                    # Initialize points in the m^th range as not xElement
                    isxElement = 0 
                    for j in range(len(xElementCombIndex)):
                        if (m == xElementCombIndex[j]):
                            isxElement = xElementCombAbund[j] 
                            # Any value not equal 0 is true. Points in the 
                            # m^th range all contain isxElement number of 
                            # xElement atoms
                    # Iterate through all points in the m^th range and decide 
                    # Whether the point is in the i^th triangle    
                    for j in range(len(gridToPoint[xVoxel][yVoxel][zVoxel][m])):
                        p = [pointX[m][gridToPoint[xVoxel][yVoxel][zVoxel][m][j]],\
    pointY[m][gridToPoint[xVoxel][yVoxel][zVoxel][m][j]], \
    pointZ[m][gridToPoint[xVoxel][yVoxel][zVoxel][m][j]]] 
                        v0p = subt(p, v[0])
                        nsq = dotProd(n, n)
                        dis = dotProd(v0p, n) / nsq**0.5 # Distance of p from the 
                                                         # plane of the triangle
                        # Decide whether the j^th point in the m^th range is 
                        # close to the i^th triangle
                        # If the distance of j^th point in the m^th range from the 
                        # plane of the i^th triangle is within the range of 
                        # [0, abs(dis)]
                        if xDis >= abs(dis):
                            a = dotProd(crossProd(v0p, v02), n) / nsq
                            if a<0:
                                continue
                            b = - dotProd(crossProd(v0p, v01), n) / nsq
                            if b<0 or a+b>1:
                                continue
                            # Then the projection of j^th point in the m^th range 
                            # on the plane of the i^th triangle is located 
                            # inside the i^th triangle
                            if isxElement:
                                triangleValueOfxElement[-1]+=isxElement
                                pointOfxElement.append(p)
                            else:
                                pointOfOtherElement.append(p)
                            for k in range(len(elementCombAbund[m])):
                                triangleValueOfAllElement[-1] += \
    elementCombAbund[m][k]
    triangleValueOfxElementNew.append(triangleValueOfxElement[-1])
    triangleValueOfAllElementNew.append(triangleValueOfAllElement[-1])
    if triangleValueOfAllElement[-1] == 0:
        triangleConcentrationOfxElement.append(0)
        triangleConcentrationOfxElementNew.append(0)
    else:
        triangleConcentrationOfxElement.append(float(triangleValueOfxElement[-1])\
                / triangleValueOfAllElement[-1])
        triangleConcentrationOfxElementNew.append(float(triangleValueOfxElement[-1])\
                / triangleValueOfAllElement[-1])
time5 = time.time()
print("Calculation is finished at " + str(time5 - time0))      

# Testing
# Test if all points are well distributed into cubic grid, whether too much 
# collision happens
#maxN = 0
#for i in range(voxelNumX):
#    for j in range(voxelNumY):
#        for k in range(voxelNumZ):
#            if (len(gridToPoint[i][j][k]) > maxN):
#                maxN = len(gridToPoint[i][j][k])
#print(maxN)
#print("The numbers of voxels in the rectangular grid along x, y and z directions \
#        are respectively %d, %d, and %d"
#      % (voxelNumX, voxelNumY, voxelNumZ))
#print(gridToPoint)

# **************************************************************
# Combine bins to make a larger bin to reduce the statistical error of 
# concentration of xElement
# **************************************************************

import math
import time
time0 = time.time()
# Declare variables
edgeNumber = binSize-1
if edgeNumber == 0:
    logBinSize = -1
else:
    logBinSize = int(math.floor(math.log(edgeNumber, 2)))
closeVertices = [[set() for j in range(vertexNumber)] for i in \
                range(logBinSize+1)] # Consist of binSize number of list. 
# The j^th elememt in the i^th list is the set of all vertices that are connected 
# to the j^th vertex through at most 2^i edges.
inBinVertices = [{i} for i in range(vertexNumber)] #All the vertices at most 
                                    #binSize-1 edges away from the i^th vertex

# The vertices 2^0 edge away from a vertex are those connected to the 
# vertex to form a triangle
if not (edgeNumber == 0):
    for i in range(triangleNumber):
        for j in range(3):
            closeVertices[0][triangle[i][j]] = closeVertices[0][triangle[i][j]] |\
{triangle[i][k] for k in range(3)}

# The vertices 2^i edges away from the j^th vertex are the set of vertices at most
# 2^(i-1) edges away from those at most 2^(i-1) edges away from the j^th vertex
for i in range(1, logBinSize+1):
    for j in range(vertexNumber):
        for k in closeVertices[i-1][j]:
            closeVertices[i][j] = closeVertices[i][j] | closeVertices[i-1][k]
time1 = time.time()
print("Vertices 2^n edges away are found at " + str(time1 - time0))
# Write edgeNumber as a binary number, whose i^th digit from the right (counting
# from 0) determines if we use closeVertices[i]
for i in range(0, logBinSize+1):
    if (edgeNumber>>i) % 2: # If the i^th digit from the right is 1, then 
                            # we will use closeVertices[i]
        for j in range(vertexNumber):
            tempVertices = set()
            for k in inBinVertices[j]:
                tempVertices = tempVertices | closeVertices[i][k]
            inBinVertices[j] = tempVertices
        #testing
        timex = time.time()
        print("Vertices {} edges away are found at {}".format(edgeNumber - \
                    ((edgeNumber>>(i+1))<<(i+1)), timex - time0))

time2 = time.time()
print("Calculation finished at " + str(time2 - time0))
if not (edgeNumber == 0):
    print("Vertex 0 is directly connected to" + str(closeVertices[0][0]))
print("Vertex 0 is connected through {} edges to {}".format(edgeNumber, \
            inBinVertices[0]))
#for i in inBinVertices[0]:
#    print("Vertex {} is directly connected to {}".format(i, closeVertices[0][i]))

# *************************************************************************
# Plot the concentration of atoms of the element of interest on the surface
# *************************************************************************
# Interpolate over the surface using point data to make color smooth
# Combine triangles to create a larger bin to reduce statistical error
vertexScalarxElementCount = [0 for i in range(len(vertexX))] 
vertexScalarAllElementCount = [0 for i in range(len(vertexX))] 
vertexScalar = [0 for i in range(len(vertexX))] 
for i in range(len(triangle)):
    triangleInBinVertices = set()
    for j in range(3):
        triangleInBinVertices = triangleInBinVertices.union(\
                inBinVertices[triangle[i][j]])
    for j in triangleInBinVertices:
        vertexScalarAllElementCount[j] = vertexScalarAllElementCount[j] \
                                         + triangleValueOfAllElement[i]
        vertexScalarxElementCount[j]= vertexScalarxElementCount[j] + \
                                      triangleValueOfxElement[i]
for i in range(len(vertexX)):
    if vertexScalarAllElementCount[i] == 0:
        vertexScalar[i] = 0
    else:
        vertexScalar[i] = float(vertexScalarxElementCount[i])/\
                          vertexScalarAllElementCount[i]*100 
                          #multiply by 100 to get the percentage

# *************************************************************************
# Visualize the concentration of atoms of the element of interest on 
# the surface in 3D animation
# *************************************************************************

# Require installing mayavi
from mayavi import mlab 
surf = mlab.pipeline.triangular_mesh_source(vertexX, vertexY, vertexZ, \
        triangleNew, scalars=vertexScalar)
surf.data.cell_data.add_array([triangleConcentrationOfxElementNew[i]*100 \
        for i in range(len(triangleConcentrationOfxElementNew))])
surf.data.cell_data.get_array(0).name = 'value of x element'
surf.data.point_data.add_array(vertexScalar)
surf.data.point_data.get_array(0).name = 'value of x element'
surf.update()
surfValueOfxElement = mlab.pipeline.set_active_attribute(surf, \
        point_scalars='value of x element')
#print(surfValueOfxElement.parent.parent.parent.module_manager)

# Display atoms near the interface 
mlab.points3d([elem[0] for elem in pointOfxElement], [elem[1] for elem in pointOfxElement], 
              [elem[2] for elem in pointOfxElement], color=(1,0,0), scale_mode='none', scale_factor=0.005)
mlab.points3d([elem[0] for elem in pointOfOtherElement], [elem[1] for elem in pointOfOtherElement],
              [elem[2] for elem in pointOfOtherElement], color=(1,1,1), scale_mode='none', scale_factor=0.005)

@mlab.animate(delay=10)
def anim():
    f = mlab.gcf()
    while 1:
        f.scene.camera.azimuth(1)
        f.scene.render()
        yield

a = anim()
mlab.pipeline.surface(surfValueOfxElement)
mlab.colorbar(title='concentration of ' + str(xElementName) + ' atoms in percent', orientation='horizontal', nb_labels=3)
mlab.show() 

    
    
    

