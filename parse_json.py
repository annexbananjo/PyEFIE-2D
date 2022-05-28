# Util functions for pyefie2d
from excitation import Planewave
import numpy as np
import json
from object import Object
from point2D import Point2D
from primitive_objects import Rectangle, Circle
from material import Material
from excitation import Planewave
from frequency import Frequency


def factory_rectangle(dictObj):
    """ Create a rectangle object from parsed json file
    """
    idObj = dictObj['id']
    pMin = Point2D(dictObj['pmin'][0], dictObj['pmin'][1])
    pMax = Point2D(dictObj['pmax'][0], dictObj['pmax'][1])
    rectangle = Rectangle(idObj, pMin, pMax)
    return rectangle
#--------------------------------------------

def factory_circle(dictObj):
    """ Create a circle object from parsed json file
    """
    idObj = dictObj['id']
    center = Point2D(dictObj['center'][0], dictObj['center'][1])
    radius = dictObj['radius']
    circObj = Circle(idObj, radius, center)
    return circObj
#--------------------------------------------

def factory_material(dictObj):
    """ Create a material object from parsed json file
    """
    idObj = dictObj['id']
    epsr = dictObj['epsr']
    sigma = dictObj['sigma']
    material = Material(idObj, epsr, sigma)
    return material
#--------------------------------------------

def factory_frequency(dictObj):
    """ Create a frequency object from parsed json file
    """
    idObj = dictObj['id']
    val = dictObj['val']
    freq = Frequency(idObj, val)
    return freq
#--------------------------------------------

def factory_planewave(dicObj):
    numPlws = dicObj['num_plws']
    startPhi = dicObj['start_phi']
    endPhi = dicObj['end_phi']
    plw = Planewave(numPlws, startPhi, endPhi)
    return plw
#--------------------------------------------

def create_object_to_material_map(fileJSON):
    """ Given a json file, it creates a map between the object and the material
        (as a dictionary)
    """
    objToMatMap = {}
    # Open the JSON file
    with open(fileJSON, 'r') as json_file:
        data = json.load(json_file)
        # loop over all the geometry entities
        for d in data['geometry']:
            if d == 'circle' or d == 'rectangle':
                # for obj in data['geometry'][d]:
                #     objToMatMap.update({obj['id']: obj['id_mat']})
                obj =data['geometry'][d]
                objToMatMap.update({obj['id']: obj['id_mat']})
        return objToMatMap
#--------------------------------------------

def create_material_to_frequency_map(fileJSON):
    """ Given a json file, it creates a map between the material and frequency
        (as a dictionary)
    """
    matToFreqMap = {}
    # Open the JSON file
    with open(fileJSON, 'r') as json_file:
        data = json.load(json_file)
        # loop over all the geometry entities
        # for obj in data['material']:
        #     matToFreqMap.update({obj['id']: obj['id_freq']})
        obj = data['material']
        matToFreqMap.update({obj['id']: obj['id_freq']})
        return matToFreqMap
#--------------------------------------------

def parse_json(fileJSON):
    """ General function that parses the json file and creates all the maps
    """
    # First loop over the geometry
    vecObj = []
    vecFreq = []
    vecMat = []
    vecInc = []
    # Open the JSON file
    with open(fileJSON, 'r') as json_file:
        data = json.load(json_file)
        # Resolution
        resolution = data['solver']['resolution']
        # print(resolution)
        # Geometry
        for d in data['geometry']:
            print(d)
            if d == 'circle':
                # for obj in data['geometry'][d]:
                #     circ = factory_circle(obj)
                #     vecObj.append(circ)
                obj = data['geometry'][d]
                circ = factory_circle(obj)
                vecObj.append(circ)
            elif d == 'rectangle':
                # for obj in data['geometry'][d]:
                #     rec = factory_rectangle(obj)
                #     vecObj.append(rec)
                obj = data['geometry'][d]
                rec = factory_rectangle(obj)
                vecObj.append(rec)
            else:
                raise ValueError("Wrong geometry entity")
        # Frequency
        # for obj in data['solver']['frequency']:
        #     freq = factory_frequency(obj)
        #     vecFreq.append(freq)
        obj =data['solver']['frequency']
        freq = factory_frequency(obj)
        vecFreq.append(freq)
        # Material
        # for obj in data['material']:
        #     mat = factory_material(obj)
        #     vecMat.append(mat)
        obj = data['material']
        mat = factory_material(obj)
        vecMat.append(mat)
        # Excitation
        for d in data['excitation']:
            if d == 'planewave':
                # for obj in data['excitation'][d]:
                #     plw = factory_planewave(obj)
                #     vecInc.append(plw)
                obj = data['excitation'][d]
                plw = factory_planewave(obj)
                vecInc.append(plw)
            else:
                raise ValueError("Wrong excitation")

    # Create the maps
    objToMatMap = create_object_to_material_map(fileJSON)
    matToFreqMap = create_material_to_frequency_map(fileJSON)
    # Return all the info
    return resolution, vecObj, vecInc, vecFreq, vecMat, objToMatMap, matToFreqMap
#--------------------------------------------


if __name__ == '__main__':
    fileJSON = r'C:\Users\mgban\Desktop\data\PetersonCylinder\PetersonCyl.json'
    # fileJSON = 'data/PetersonCylinder/PetersonCyl.json'
    resolution, vecObj, vecInc, vecFreq, vecMat, objToMatMap, matToFreqMap = parse_json(
        fileJSON)
