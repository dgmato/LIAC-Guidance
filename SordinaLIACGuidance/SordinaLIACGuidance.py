import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

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

    logic = SordinaLIACGuidanceLogic()

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
    self.upDownRecordButton.connect('clicked(bool)', self.onUpDownRecordButton)
    self.upDownStopButton.connect('clicked(bool)', self.onUpDownStopButton)
    self.upDownRestartButton.connect('clicked(bool)', self.onUpDownRestartButton)
    self.rollRecordButton.connect('clicked(bool)', self.onRollRecordButton)
    self.rollStopButton.connect('clicked(bool)', self.onRollStopButton)
    self.rollRestartButton.connect('clicked(bool)', self.onRollRestartButton)
    
    
    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    pass

  def onUpDownRecordButton(self):
    self.upDownRecordButton.enabled = False
    self.upDownStopButton.enabled = True
    self.upDownRestartButton.enabled = False
    

  def onUpDownStopButton(self):
    self.upDownRecordButton.enabled = False
    self.upDownStopButton.enabled = False
    self.upDownRestartButton.enabled = True
    

  def onUpDownRestartButton(self):
    self.upDownRecordButton.enabled = True
    self.upDownStopButton.enabled = False
    self.upDownRestartButton.enabled = False
    

  def onRollRecordButton(self):
    self.rollRecordButton.enabled = False
    self.rollStopButton.enabled = True
    self.rollRestartButton.enabled = False
    
  def onRollStopButton(self):
    self.rollRecordButton.enabled = False
    self.rollStopButton.enabled = False
    self.rollRestartButton.enabled = True

  def onRollRestartButton(self):
    self.rollRecordButton.enabled = True
    self.rollStopButton.enabled = False
    self.rollRestartButton.enabled = False

#
# SordinaLIACGuidanceLogic
#

class SordinaLIACGuidanceLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() == None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : imageThreshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('SordinaLIACGuidanceTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


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

  def test_SordinaLIACGuidance1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = SordinaLIACGuidanceLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
