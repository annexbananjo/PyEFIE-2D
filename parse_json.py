# Util functions for pyefie2d
import numpy as np
import json
from object import Object
from point2D import Point2D
from primitive_objects import Rectangle, Circle
from material import Material
from frequency import Frequency


def factory_rectangle(dictObj):
    """ Create a rectangle object from parsed json file
    """
    idObj = dictObj['id']
    pMin = Point2D(dictObj['pmin'][0], dictObj['pmin'][1])
    pMax = Point2D(dictObj['pmax'][0], dictObj['pmax'][1])
    rectangle = Rectangle(idObj, pMin, pMax)
    return rectangle


def factory_circle(dictObj):
    """ Create a circle object from parsed json file
    """
    idObj = dictObj['id']
    center = Point2D(dictObj['center'][0], dictObj['center'][1])
    radius = dictObj['radius']
    circObj = Circle(idObj, radius, center)
    return circObj


def factory_material(dictObj):
    """ Create a material object from parsed json file
    """
    idObj = dictObj['id']
    epsr = dictObj['epsr']
    sigma = dictObj['sigma']
    material = Material(idObj, epsr, sigma)
    return material


def factory_frequency(dictObj):
    """ Create a frequency object from parsed json file
    """
    idObj = dictObj['id']
    val = dictObj['val']
    freq = Frequency(idObj, val)
    return freq


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
                for obj in data['geometry'][d]:
                    objToMatMap.update({obj['id']: obj['id_mat']})
        return objToMatMap


def create_material_to_frequency_map(fileJSON):
    """ Given a json file, it creates a map between the material and frequency
        (as a dictionary)
    """
    matToFreqMap = {}
    # Open the JSON file
    with open(fileJSON, 'r') as json_file:
        data = json.load(json_file)
        # loop over all the geometry entities
        for obj in data['material']:
            matToFreqMap.update({obj['id']: obj['id_freq']})
        return matToFreqMap


def parse_json(fileJSON):
    """ General function that parses the json file and creates all the maps 
    """
    # First loop over the geometry
    vecObj = []
    vecFreq = []
    vecMat = []
    # Open the JSON file
    with open(fileJSON, 'r') as json_file:
        data = json.load(json_file)
        # Resolution
        resolution = data['resolution']
        # Geometry
        for d in data['geometry']:
            if d == 'circle':
                for obj in data['geometry'][d]:
                    circ = factory_circle(obj)
                    vecObj.append(circ)
            elif d == 'rectangle':
                for obj in data['geometry'][d]:
                    rec = factory_rectangle(obj)
                    vecObj.append(rec)
            else:
                raise ValueError("Wrong geometry entity")
        # Frequency
        for obj in data['frequency']:
            freq = factory_frequency(obj)
            vecFreq.append(freq)
        # Material
        for obj in data['material']:
            mat = factory_material(obj)
            vecMat.append(mat)
        # Excitation is missing!
    # Create the maps
    objToMatMap = create_object_to_material_map(fileJSON)
    matToFreqMap = create_material_to_frequency_map(fileJSON)
    # Return all the info
    return resolution, vecObj, vecFreq, vecMat, objToMatMap, matToFreqMap


if __name__ == '__main__':
    # fileJSON = 'data/ThreeObjects/ThreeObjects.json'
    fileJSON = 'data/PetersonCylinder/PetersonCyl.json'
    resolution, vecObj, vecFreq, vecMat, objToMatMap, matToFreqMap = parse_json(
        fileJSON)
    a = 1
