import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy
import time, csv

# import SordinaLIAC
# import PyRobotics


#
# SordinaLIACGuidance
#

class SordinaLIACGuidance(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "SordinaLIACGuidance" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Eugenio Marinetto (LIM)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# SordinaLIACGuidanceWidget
#

class SordinaLIACGuidanceWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    ########################################################
    ############# Transform Definition Area ################
    ########################################################
    transformsCollapsibleButton = ctk.ctkCollapsibleButton()
    transformsCollapsibleButton.text = "Transform Definition"
    self.layout.addWidget(transformsCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(transformsCollapsibleButton)

    # ApplicatorToTracker transform selector
    self.applicatorToTrackerSelector = slicer.qMRMLNodeComboBox()
    self.applicatorToTrackerSelector.nodeTypes = ( ("vtkMRMLLinearTransformNode"), "" )
    self.applicatorToTrackerSelector.selectNodeUponCreation = True
    self.applicatorToTrackerSelector.addEnabled = False
    self.applicatorToTrackerSelector.removeEnabled = False
    self.applicatorToTrackerSelector.noneEnabled = False
    self.applicatorToTrackerSelector.showHidden = False
    self.applicatorToTrackerSelector.showChildNodeTypes = False
    self.applicatorToTrackerSelector.setMRMLScene( slicer.mrmlScene )
    self.applicatorToTrackerSelector.setToolTip( "Pick the applicatorToTracker transform." )
    parametersFormLayout.addRow("ApplicatorToTracker transform: ", self.applicatorToTrackerSelector)
 
    # LiacToTracker transform selector
    self.liacToTrackerSelector = slicer.qMRMLNodeComboBox()
    self.liacToTrackerSelector.nodeTypes = ( ("vtkMRMLLinearTransformNode"), "" )
    self.liacToTrackerSelector.selectNodeUponCreation = True
    self.liacToTrackerSelector.addEnabled = False
    self.liacToTrackerSelector.removeEnabled = False
    self.liacToTrackerSelector.noneEnabled = False
    self.liacToTrackerSelector.showHidden = False
    self.liacToTrackerSelector.showChildNodeTypes = False
    self.liacToTrackerSelector.setMRMLScene( slicer.mrmlScene )
    self.liacToTrackerSelector.setToolTip( "Pick the LiacToTracker transform." )
    parametersFormLayout.addRow("LiacToTracker transform: ", self.liacToTrackerSelector)

    ########################################################
    ################## Calibration Area ####################
    ########################################################
    calibrationCollapsibleButton = ctk.ctkCollapsibleButton()
    calibrationCollapsibleButton.text = "Calibration"
    self.layout.addWidget(calibrationCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(calibrationCollapsibleButton)

    #
    # Icons retrieval
    #
    sordinaLIACGuidanceModuleDirectoryPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","")
    recordIcon = qt.QIcon(sordinaLIACGuidanceModuleDirectoryPath + '/Resources/Icons/recordIcon.png')
    stopIcon = qt.QIcon(sordinaLIACGuidanceModuleDirectoryPath + '/Resources/Icons/stopIcon.png')
    restartIcon = qt.QIcon(sordinaLIACGuidanceModuleDirectoryPath + '/Resources/Icons/restartIcon.png')
    saveIcon = qt.QIcon(sordinaLIACGuidanceModuleDirectoryPath + '/Resources/Icons/saveIcon.png')    

    #
    # Up-down movement recording
    #
    upDownMovementGroupBox = ctk.ctkCollapsibleGroupBox()
    upDownMovementGroupBox.setTitle("Up-Down Movement")
    upDownMovementGroupBoxFormLayout = qt.QFormLayout(upDownMovementGroupBox)
    parametersFormLayout.addRow(upDownMovementGroupBox)
  
    upDownMovementLayout = qt.QHBoxLayout()
    upDownMovementGroupBoxFormLayout.addRow(upDownMovementLayout)

    self.upDownRecordButton = qt.QPushButton(" Record")
    self.upDownRecordButton.setIcon(recordIcon)
    self.upDownRecordButton.enabled = True
    upDownMovementLayout.addWidget(self.upDownRecordButton)
    # upDownMovementGroupBoxFormLayout.addRow(self.upDownRecordButton)

    self.upDownStopButton = qt.QPushButton(" Stop")
    self.upDownStopButton.setIcon(stopIcon)
    self.upDownStopButton.enabled = False
    upDownMovementLayout.addWidget(self.upDownStopButton)
    # upDownMovementGroupBoxFormLayout.addRow(self.upDownStopButton)

    self.upDownRestartButton = qt.QPushButton(" Restart")
    self.upDownRestartButton.setIcon(restartIcon)
    self.upDownRestartButton.enabled = False
    upDownMovementLayout.addWidget(self.upDownRestartButton)
    # upDownMovementGroupBoxFormLayout.addRow(self.upDownRestartButton)

    #
    # Roll movement recording
    #
    rollMovementGroupBox = ctk.ctkCollapsibleGroupBox()
    rollMovementGroupBox.setTitle("Roll Movement")
    rollMovementGroupBoxFormLayout = qt.QFormLayout(rollMovementGroupBox)
    parametersFormLayout.addRow(rollMovementGroupBox)

    rollMovementLayout = qt.QHBoxLayout()
    rollMovementGroupBoxFormLayout.addRow(rollMovementLayout)
  
    self.rollRecordButton = qt.QPushButton(" Record")
    self.rollRecordButton.setIcon(recordIcon)
    self.rollRecordButton.enabled = True
    rollMovementLayout.addWidget(self.rollRecordButton)
    # rollMovementGroupBoxFormLayout.addRow(self.rollRecordButton)

    self.rollStopButton = qt.QPushButton(" Stop")
    self.rollStopButton.setIcon(stopIcon)
    self.rollStopButton.enabled = False
    rollMovementLayout.addWidget(self.rollStopButton)
    # rollMovementGroupBoxFormLayout.addRow(self.rollStopButton)

    self.rollRestartButton = qt.QPushButton(" Restart")
    self.rollRestartButton.setIcon(restartIcon)
    self.rollRestartButton.enabled = False
    rollMovementLayout.addWidget(self.rollRestartButton)
    # rollMovementGroupBoxFormLayout.addRow(self.rollRestartButton)

    #
    # Save/load calibration file
    #
    self.saveCalibrationButton = qt.QPushButton(" Save Calibration")
    self.saveCalibrationButton.setIcon(saveIcon)
    self.saveCalibrationButton.enabled = True
    self.saveCalibrationButton.adjustSize()
    
    self.loadCalibrationButton = qt.QPushButton(" Load Calibration")
    self.loadCalibrationButton.setIcon(saveIcon)
    self.loadCalibrationButton.enabled = True
    self.loadCalibrationButton.adjustSize()

    self.calibrationErrorLabel = qt.QLabel('Calibration error: - mm')
    self.calibrationErrorLabel.setStyleSheet("QLabel { color : #000000; font: bold}" )
    
    parametersFormLayout.addRow(self.loadCalibrationButton, self.saveCalibrationButton)
    parametersFormLayout.addRow(self.calibrationErrorLabel)  
    
    #########################################################
    #################### Guidance Area ######################
    #########################################################
    guidanceCollapsibleButton = ctk.ctkCollapsibleButton()
    guidanceCollapsibleButton.text = "Guidance"
    self.layout.addWidget(guidanceCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(guidanceCollapsibleButton)

    # Load Models Button
    self.loadModelsButton = qt.QPushButton()
    self.loadModelsButton.toolTip = "Soft tissue visibility."
    self.loadModelsButton.enabled = True
    self.loadModelsButton.text = "Load 3D Models"
    parametersFormLayout.addRow(self.loadModelsButton)

    # connections
    self.applicatorToTrackerSelector.connect('currentNodeChanged(vtkMRMLNode*)',self.onSelect)
    self.liacToTrackerSelector.connect('currentNodeChanged(vtkMRMLNode*)',self.onSelect)
    self.upDownRecordButton.connect('clicked(bool)', self.onUpDownRecordButton)
    self.upDownStopButton.connect('clicked(bool)', self.onUpDownStopButton)
    self.upDownRestartButton.connect('clicked(bool)', self.onUpDownRestartButton)
    self.rollRecordButton.connect('clicked(bool)', self.onRollRecordButton)
    self.rollStopButton.connect('clicked(bool)', self.onRollStopButton)
    self.rollRestartButton.connect('clicked(bool)', self.onRollRestartButton)
    self.loadModelsButton.connect('clicked(bool)', self.onLoadModelsButton)
    
    self.SordinaLIACGuidanceLogic = SordinaLIACGuidanceLogic()  

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applicatorToTrackerTransform=self.applicatorToTrackerSelector.currentNode()
    self.liacToTrackerTransform=self.liacToTrackerSelector.currentNode()

  def onUpDownRecordButton(self):
    self.upDownRecordButton.enabled = False
    self.upDownStopButton.enabled = True
    self.upDownRestartButton.enabled = False
    
    if not self.SordinaLIACGuidanceLogic.observedNode:
      self.SordinaLIACGuidanceLogic.resetScene()
      self.SordinaLIACGuidanceLogic.record=True
      self.SordinaLIACGuidanceLogic.movementType='Up-Down'
      self.SordinaLIACGuidanceLogic.removeUpdateObserver()
      self.SordinaLIACGuidanceLogic.addUpdateObserver(self.liacToTrackerTransform)

  def onUpDownStopButton(self):
    self.upDownRecordButton.enabled = False
    self.upDownStopButton.enabled = False
    self.upDownRestartButton.enabled = True

    self.SordinaLIACGuidanceLogic.record=False
    self.SordinaLIACGuidanceLogic.saveData('Up-Down')
    self.SordinaLIACGuidanceLogic.observedNode=None
    
  def onUpDownRestartButton(self):
    self.upDownRecordButton.enabled = True
    self.upDownStopButton.enabled = False
    self.upDownRestartButton.enabled = False
    
    self.SordinaLIACGuidanceLogic.resetScene()

  def onRollRecordButton(self):
    self.rollRecordButton.enabled = False
    self.rollStopButton.enabled = True
    self.rollRestartButton.enabled = False

    if not self.SordinaLIACGuidanceLogic.observedNode:
      self.SordinaLIACGuidanceLogic.resetScene()
      self.SordinaLIACGuidanceLogic.record=True
      self.SordinaLIACGuidanceLogic.movementType='Roll'
      self.SordinaLIACGuidanceLogic.removeUpdateObserver()
      self.SordinaLIACGuidanceLogic.addUpdateObserver(self.liacToTrackerTransform)
    
  def onRollStopButton(self):
    self.rollRecordButton.enabled = False
    self.rollStopButton.enabled = False
    self.rollRestartButton.enabled = True

    self.SordinaLIACGuidanceLogic.record=False
    self.SordinaLIACGuidanceLogic.saveData('Roll')
    self.SordinaLIACGuidanceLogic.observedNode=None

  def onRollRestartButton(self):
    self.rollRecordButton.enabled = True
    self.rollStopButton.enabled = False
    self.rollRestartButton.enabled = False

    self.SordinaLIACGuidanceLogic.resetScene()

  def onLoadModelsButton(self):
    # Load head model + Apply transform
    self.applicator1Model = slicer.util.getNode('Applicator1_origin')
    if not self.applicator1Model:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'Applicator1_origin.stl')
        self.applicator1Model = slicer.util.getNode(pattern="Applicator1_origin")
        self.applicator1ModelDisplay=self.applicator1Model.GetModelDisplayNode()
        self.applicator1ModelDisplay.SetColor([1,0.7,0.53])

    self.applicator2Model = slicer.util.getNode('Applicator2_origin')
    if not self.applicator2Model:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'Applicator2_origin.stl')
        self.applicator2Model = slicer.util.getNode(pattern="Applicator2_origin")
        self.applicator2ModelDisplay=self.applicator2Model.GetModelDisplayNode()
        self.applicator2ModelDisplay.SetColor([1,0.7,0.53])

    self.armRotateModel = slicer.util.getNode('ArmRotate_origin')
    if not self.armRotateModel:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'ArmRotate_origin.stl')
        self.armRotateModel = slicer.util.getNode(pattern="ArmRotate_origin")
        self.armRotateModelDisplay=self.armRotateModel.GetModelDisplayNode()
        self.armRotateModelDisplay.SetColor([1,0.7,0.53])

    self.armRotateFrontModel = slicer.util.getNode('ArmRotateFront_origin')
    if not self.armRotateFrontModel:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'ArmRotateFront_origin.stl')
        self.armRotateFrontModel = slicer.util.getNode(pattern="ArmRotateFront_origin")
        self.armRotateFrontModelDisplay=self.armRotateFrontModel.GetModelDisplayNode()
        self.armRotateFrontModelDisplay.SetColor([1,0.7,0.53])

    self.bodyModel = slicer.util.getNode('Body_origin')
    if not self.bodyModel:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'Body_origin.stl')
        self.bodyModel = slicer.util.getNode(pattern="Body_origin")
        self.bodyModelDisplay=self.bodyModel.GetModelDisplayNode()
        self.bodyModelDisplay.SetColor([1,0.7,0.53])

    self.elevationModel = slicer.util.getNode('Elevation_origin')
    if not self.elevationModel:
        sordinaLIACGuidanceModuleDataPath = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'Resources/Models/'
        slicer.util.loadModel(sordinaLIACGuidanceModuleDataPath + 'Elevation_origin.stl')
        self.elevationModel = slicer.util.getNode(pattern="Elevation_origin")
        self.elevationModelDisplay=self.elevationModel.GetModelDisplayNode()
        self.elevationModelDisplay.SetColor([1,0.7,0.53])

#
# SordinaLIACGuidanceLogic
#

class SordinaLIACGuidanceLogic(ScriptedLoadableModuleLogic):
  
  def __init__(self): 
    self.m = vtk.vtkMatrix4x4() # 4x4 VTK matrix to save the transformation matrix sent through Plus.
    self.observedNode = None
    self.outputObserverTag = -1
    self.record = False
    self.movementType = ''
    self.timerActive=False
    self.myTimer=Timer()
    self.recordedDataBuffer = [] 
    self.timeStamp = numpy.array([])
    self.positionX = numpy.array([])
    self.positionY = numpy.array([])
    self.positionZ = numpy.array([])
    
  def resetScene(self):
    self.m = vtk.vtkMatrix4x4() # 4x4 VTK matrix to save the transformation matrix sent through Plus.
    self.observedNode = None
    self.outputObserverTag = -1
    self.record = False
    self.movementType = ''
    self.timerActive=False
    self.myTimer=Timer()
    self.recordedDataBuffer = [] 
    self.timeStamp = numpy.array([])
    self.positionX = numpy.array([])
    self.positionY = numpy.array([])
    self.positionZ = numpy.array([])
      
  def addUpdateObserver(self, inputNode):
    """
    Summary: Add update observer to inputNode.
    """
    self.observedNode = inputNode
    if self.outputObserverTag == -1:
      self.outputObserverTag = inputNode.AddObserver('ModifiedEvent', self.updateSceneCallback)

  def removeUpdateObserver(self):
    """
    Summary: Remove update observer.
    """
    if self.outputObserverTag != -1:
      self.observedNode.RemoveObserver(self.outputObserverTag)
      self.outputObserverTag = -1
      self.observedNode = None

  def updateSceneCallback(self, modifiedNode, event=None): 
    """
    Summary: This functions is called when the observed node (to which an observer has been added) is modified.
    """           
    if self.record:
        if self.timerActive==False:
            self.myTimer.startTimer()
            self.timerActive=True
        self.recordPositionsForCalibration(self.observedNode)

  def recordPositionsForCalibration(self, InputTransformMatrix):
    print('Recording...')
    # Input tranformation node is saved into a VTK matrix 'm'
    InputTransformMatrix.GetMatrixTransformToParent(self.m)
    # VTK matrix containing the input data is tranformed into a ordinary numpy matrix
    self.positionX = numpy.append(self.positionX, self.m.GetElement(0, 3))
    self.positionY = numpy.append(self.positionY, self.m.GetElement(1, 3))
    self.positionZ = numpy.append(self.positionZ, self.m.GetElement(2, 3))
    
    t=self.myTimer.getElapsedTime()
    self.timeStamp=numpy.append(self.timeStamp, t) # Time stamp saved for each iteration.
  
  def saveData(self, movementType):
    print('Saving data into .csv file...')
    self.myTimer.stopTimer()
    dateAndTime = time.strftime("_%Y-%m-%d_%H-%M-%S")
    if movementType=='Up-Down':
      path = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'SavedData/' + 'CalibrationData_UpDown_' + dateAndTime
    elif movementType=='Roll':
      path = slicer.modules.sordinaliacguidance.path.replace("SordinaLIACGuidance.py","") + 'SavedData/' + 'CalibrationData_Roll_' + dateAndTime
    else:
      print ('ERROR saveData method: non-valid movement type.')
    # Save data to .csv 
    self.writeRecordedDataBufferToFile(path + '.csv')  
    
  def writeRecordedDataBufferToFile(self, path): 
    with open(path, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=",")
      writer.writerow(['timeStamp', 'x' , 'y', 'z'])

      for i in range(len(self.timeStamp)):
        vector = numpy.array([self.timeStamp[i],\
            self.positionX[i],\
            self.positionY[i],\
            self.positionZ[i]])
        writer.writerow(vector)


class SordinaLIACGuidanceTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_SordinaLIACGuidance1()


class Timer(object):

  def __init__(self):
    self.startTime = 0.0
    self.stopTime = 0.0
    self.timerStarted = False
    
  def startTimer(self):
    if not self.timerStarted:      
      self.startTime = time.clock()
      if self.stopTime != 0.0:
        self.stopTime = time.clock() - self.stopTime
      self.timerStarted = True
    else:
      logging.warning('Timer already running')
      
  def stopTimer(self):
    if self.timerStarted:
      self.stopTime = time.clock()
      self.timerStarted = False
    else:
      logging.warning('Timer not running')
      
  def getElapsedTime(self):
    if self.startTime == 0.0:
      return 0.0
    elif self.stopTime == 0.0:
      return time.clock() - self.startTime
    else:
      return time.clock() - (self.stopTime - self.startTime)
        
  def resetTimer(self):
    if self.startTime != 0.0:
      self.startTime = time.clock()    
      self.stopTime = 0.0