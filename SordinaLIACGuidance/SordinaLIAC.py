import PyRobotics
import copy
import numpy
import scipy
import plotly
import vtk


class LIAC(PyRobotics.RobotArm):
    """@brief Class representing a chain of joints of the LIAC Sordina

    Initialization is made from Denavit-Hartenberg parameters of the joint inside a chain
    """

    def __init__(self, dApp = 0.5):
        """
            ## Denavit Hartenberg parameters

            | Art | Rot Z | Tra Z     | Tra X | Rot X | Min q | Max q |
            |-----|-------|-----------|-------|-------|-------|-------|
            | 1   | 0     | 0         | 0     | 0     | -90   | 90    |
            | 2   | 90    | 0         | 0     | 90    | 0     | 0.5   |
            | 3   | 0     | 0         | 0     | -90   | -90   | 90    |
            | 4   | 0     | 0         | 0     | 90    | 0     | 0.5   |
            | 5   | 0     | 0         | 0     | -90   | 0     | 0.89  |
            | 6   | 0     | d3        | 0     | 90    | -60   | 60    |
            | 7   | -90   | d4        | 0     | -90   | -30   | 15    |
            | 8   | -90   | 0         | d5    | -90   | -     | -     | Fake to reach the end point
            | 9   | 0     | d6 + dApp | 0     | 180   | -     | -     | Fake to reach the end point

            @brief create a robotArm
        """

        PyRobotics.RobotArm.__init__(self)

        # In metres
        self.d1 = 0.760/2.0
        self.d2 = 1.200
        self.d3 = 1.325
        self.d4 = 0.300
        self.d5 = 0.455
        self.d6 = 0.430

        ## Joint definition
        self.Joint1 = PyRobotics.Joint(   0,   0,  0,   0, 'rotation', idx = 1) # A, alpha, D, offset, freedom, idx
        self.Joint2 = PyRobotics.Joint(   0,  90,  0,  90, 'periscopic', idx = 2)
        self.Joint3 = PyRobotics.Joint(   0, -90,  0,   0, 'rotation', idx = 3)
        self.Joint4 = PyRobotics.Joint(   0,  90,  0,   0, 'periscopic', idx = 4)
        self.Joint5 = PyRobotics.Joint(   0, -90,  0,   0, 'periscopic', idx = 5)
        self.Joint6 = PyRobotics.Joint(   0,  90, self.d3,   0, 'rotation', idx = 6)
        self.Joint7 = PyRobotics.Joint(   0, -90, self.d4, -90, 'rotation', idx = 7)
        self.Joint8 = PyRobotics.Joint(  self.d5, -90,  0, -90, 'periscopic', idx = 8)
        self.Joint9 = PyRobotics.Joint(   0,   0,  self.d6 + dApp, 180, 'rotation', idx = 9)

        self.Joint1.setLimits(-90,90)
        self.Joint2.setLimits(-1.5,1.5)
        self.Joint3.setLimits(-30,30)
        self.Joint4.setLimits(-0.5,0.5)
        self.Joint5.setLimits(0,0.89)
        self.Joint6.setLimits(-60,60)
        self.Joint7.setLimits(-30,15)
        self.Joint8.setLimits(0,0)
        self.Joint9.setLimits(-180,180) # Fake can rotate as need

        # Load Models
        rootPath = '/Users/nenetto/Dev/liac/LiacSordinaSuperBuild/PyRobotics/'

        self.Joint1.setModel(pathToSTLModel = '', color = [0,1,1], opacity = 0.0, length = 1.0, radius = 0.1)
        self.Joint2.setModel(pathToSTLModel = '', color = [0,1,1], opacity = 0.0, length = 1.0, radius = 0.1)
        self.Joint3.setModel(pathToSTLModel = '', color = [0,1,1], opacity = 0.0, length = 1.0, radius = 0.1)

        # Body of LIAC
        self.Joint4.setModel(pathToSTLModel = rootPath + 'Body_origin.stl', color = [0,0.5,0.5], opacity = 0.4)

        # UpDown Segment
        self.Joint5.setModel(pathToSTLModel = rootPath + 'Elevation_origin.stl', color = [1,1,0], opacity = 0.4)

        # Rotation Segment (Roll)
        self.Joint6.setModel(pathToSTLModel = rootPath + 'ArmRotate_origin.stl', color = [0,1,1], opacity = 0.4)

        # Rotation Segment Front (Pitch)
        self.Joint7.setModel(pathToSTLModel = rootPath + 'ArmRotateFront_origin.stl', color = [0.5,1,0], opacity = 0.4)

        # Applicator first END (always fixed) (Change if Applicator change the diameter)
        self.Joint8.setModel(pathToSTLModel = rootPath + 'Applicator1_origin.stl', color = [0,1,0.5], opacity = 0.4)

        # Applicator first END (always fixed) (Change if Applicator change the diameter or besel)
        self.Joint9.setModel(pathToSTLModel = rootPath + 'Applicator2_origin.stl', color = [1,0,0], opacity = 0.4)

        self.pathToModelApplicator = rootPath + 'Applicator2_origin.stl'

        self.connectJoint(self.Joint1)
        self.connectJoint(self.Joint2)
        self.connectJoint(self.Joint3)
        self.connectJoint(self.Joint4)
        self.connectJoint(self.Joint5)
        self.connectJoint(self.Joint6)
        self.connectJoint(self.Joint7)
        self.connectJoint(self.Joint8)
        self.connectJoint(self.Joint9)

        for j in self.JointList:
            j.setAxesLength(axesLength = 0.2)

        self._rot1 = 0.0
        self._tra1 = 0.0
        self._rot2 = 0.0
        self._tra2 = 0.0
        self._UpDown = 0.0
        self._Roll = 0.0
        self._Pitch = 0.0

        self.MaximumRange = None

        self.calibrated = False
        self.objetiveFixed = False

    def Rotate1(self,q):
        self._rot1 += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def Translate1(self,q):
        self._tra1 += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def Rotate2(self,q):
        self._rot2 += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def Translate2(self,q):
        self._tra2 += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def UpDown(self,q):
        self._UpDown += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def Roll(self,q):
        self._Roll += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def Pitch(self,q):
        self._Pitch += q
        self.moveRobot([self._rot1,self._tra2,self._rot2,self._tra2,self._UpDown,self._Roll,self._Pitch,0,0])

    def MoveToOrigin(self):

        self.moveRobot([0,0,0,0,0,0,0,0,0])

    def ChangeApplicator(self, applicator = 0.2):

        self.dApp = applicator

        self.Joint9 = PyRobotics.Joint(   0,   0,  self.d6 + self.dApp, 180, 'rotation', idx = 9)
        self.Joint9.setLimits(-180,180)
        self.Joint9.setFrame(self.JointList[8])
        self.JointList[9] = self.Joint9

    def calibrateOpticalToLIACend(self, pointsQ5Optical, pointsQ6Optical, dAppCalibration = 0.0):
        tcalibration = pointsQ5Optical[0,0:3]

        P = pointsQ5Optical[:,0:3].copy()
        center =  numpy.mean(P,axis=0)
        Pcentered = P
        for i in range(P.shape[0]):
            Pcentered[i,:] = Pcentered[i,:] -center

        U, s, Vt = numpy.linalg.svd(numpy.dot(Pcentered.T,Pcentered), full_matrices=False)
        Vz0 = Vt.T[:,0]

        # Check direction
        if numpy.dot(Pcentered[0,:],Vz0) < numpy.dot(Pcentered[-1,:],Vz0):
            Vz0 = -Vz0

        P = pointsQ6Optical[:,0:3].copy()
        Pcentered = P
        Pprojected = P.copy()
        for i in range(P.shape[0]):
            Pcentered[i,:] -= tcalibration
            Pprojected[i,:] = Pcentered[i,:] - numpy.dot(Pcentered[i,:],Vz0)*Vz0

        U, s, Vt = numpy.linalg.svd(numpy.dot(Pprojected.T,Pprojected), full_matrices=False)
        Vy0 = Vt.T[:,0]

        if numpy.dot(Pcentered[0,:],Vy0) > numpy.dot(Pcentered[1,:],Vy0):
            Vy0 = -Vy0


        Taxis = numpy.eye(4)
        Taxis[0:3,0] = numpy.cross(Vy0,Vz0)
        Taxis[0:3,1] = Vy0
        Taxis[0:3,2] = Vz0
        Taxis[0:3,3] = tcalibration[0:3]

        Taxis[2,3] -= dAppCalibration

        self.OpticalToLIACend = Taxis.copy()


        # Place LIAC to the resting state
        self.moveRobot([0,0,0,0,0,0,0,0,0])

        self.TWorldToOptical = numpy.dot(self.JointList[-1].Axis,numpy.linalg.inv(self.OpticalToLIACend))
        self.TOpticalToWorld = numpy.linalg.inv(self.TWorldToOptical)

        self.TJoint4ToOptical = numpy.dot(numpy.linalg.inv(self.Joint4.Axis), self.TWorldToOptical)

        ## Add the Optical Cameras Axis
        self.moveRobot([0,0,0,0,0,0,0,0,0])
        self.TOptical = vtk.vtkTransform()
        self.TOptical.SetMatrix([self.TWorldToOptical[0,0],self.TWorldToOptical[0,1],self.TWorldToOptical[0,2],self.TWorldToOptical[0,3],
                        self.TWorldToOptical[1,0],self.TWorldToOptical[1,1],self.TWorldToOptical[1,2],self.TWorldToOptical[1,3],
                        self.TWorldToOptical[2,0],self.TWorldToOptical[2,1],self.TWorldToOptical[2,2],self.TWorldToOptical[2,3],
                        self.TWorldToOptical[3,0],self.TWorldToOptical[3,1],self.TWorldToOptical[3,2],self.TWorldToOptical[3,3]])
        self.Opticalaxes = vtk.vtkAxesActor()
        self.Opticalaxes.SetUserTransform(self.TOptical)
        self.Opticalaxes.SetXAxisLabelText('X Optical')
        self.Opticalaxes.SetYAxisLabelText('Y Optical')
        self.Opticalaxes.SetZAxisLabelText('Z Optical')
        self.Opticalaxes.SetTotalLength(0.2,0.2,0.2)

        self.calibrated = True

    def computeMaximumRange(self):
        """
        @brief computeMaximumRange: Calculate the volume of robot actuation and return a drawable (plotly) object

        @return drawable robot squema + surface of the boundaries
        """

        if self.MaximumRange is None:

            # Create number of points of every joint to be simulated
            N = 4 # More than 5 is very heavy for computation

            # Create simulation vectors

            paramsQrange = numpy.zeros((self.JointNumber,N))

            for i in list([1,2,3,4,5,6,7]):
                joint = self.getJoint(i)
                paramsQrange[i,:] = numpy.linspace(joint.minParam,joint.maxParam,N)

            paramsQ = numpy.meshgrid(*paramsQrange[:-1])

            Nsimulations = len(paramsQ[0].flatten())
            rangelocations = numpy.zeros((Nsimulations,3))

            for i in range(Nsimulations):
                q = list()
                for j in range(self.JointNumber -1):
                    q.append(paramsQ[j].flatten()[i])
                q.append(0.0) # For the fake joint
                self.moveRobot(q)
                rangelocations[i,:]= self.getRobotEnd()[0:3,3]

            surfaceHull =   plotly.graph_objs.Mesh3d(
                                                        name = "Range Limits",
                                                        x = rangelocations[:,0],
                                                        y = rangelocations[:,1],
                                                        z = rangelocations[:,2],
                                                        opacity=0.3,
                                                        color = 'yellow',
                                                        alphahull = 0
                                                      )
            self.MaximumRange = surfaceHull
            self.RangeLocations = rangelocations

        return self.MaximumRange

    def SetObjective(self, ApplicatorObjetiveInOpticalCoordinates):

        if self.calibrated:
            # Save the Transformation Goal
            self.ApplicatorObjetiveAxis = ApplicatorObjetiveInOpticalCoordinates.copy()

            self.ApplicatorObjetiveAxis = numpy.dot(self.TWorldToOptical,self.ApplicatorObjetiveAxis)

            self.ApplicatorObjetive = vtk.vtkTransform()
            self.ApplicatorObjetive.SetMatrix([self.ApplicatorObjetiveAxis[0,0],self.ApplicatorObjetiveAxis[0,1],self.ApplicatorObjetiveAxis[0,2],self.ApplicatorObjetiveAxis[0,3],
                            self.ApplicatorObjetiveAxis[1,0],self.ApplicatorObjetiveAxis[1,1],self.ApplicatorObjetiveAxis[1,2],self.ApplicatorObjetiveAxis[1,3],
                            self.ApplicatorObjetiveAxis[2,0],self.ApplicatorObjetiveAxis[2,1],self.ApplicatorObjetiveAxis[2,2],self.ApplicatorObjetiveAxis[2,3],
                            self.ApplicatorObjetiveAxis[3,0],self.ApplicatorObjetiveAxis[3,1],self.ApplicatorObjetiveAxis[3,2],self.ApplicatorObjetiveAxis[3,3]])
            self.ApplicatorObjetiveaxes = vtk.vtkAxesActor()
            self.ApplicatorObjetiveaxes.SetUserTransform(self.ApplicatorObjetive)
            self.ApplicatorObjetiveaxes.SetXAxisLabelText('X App')
            self.ApplicatorObjetiveaxes.SetYAxisLabelText('Y App')
            self.ApplicatorObjetiveaxes.SetZAxisLabelText('Z App')
            self.ApplicatorObjetiveaxes.SetTotalLength(0.2,0.2,0.2)

            self.STLReaderApplicator = vtk.vtkSTLReader()
            self.STLReaderApplicator.SetFileName(self.pathToModelApplicator)
            self.STLReaderApplicator.Update()

            self.mapperApplicator = vtk.vtkPolyDataMapper()
            self.mapperApplicator.SetInputConnection(self.STLReaderApplicator.GetOutputPort())

            self.modelActorApplicator = vtk.vtkActor()
            self.modelActorApplicator.SetMapper(self.mapperApplicator)

            self.modelActorApplicator.SetUserTransform(self.ApplicatorObjetive)
            self.modelActorApplicator.GetProperty().SetColor(0,1,0)
            self.modelActorApplicator.GetProperty().SetOpacity(0.5)

            self.objetiveFixed = True

    def solveObjective(self, draw = False, animate = False, fold = False, qInit = None):

        if qInit is None:
            qInit = numpy.zeros((self.JointNumber))

        qSol, poseError, angleError =  self.inverseKinematics(endMatrix = self.ApplicatorObjetiveAxis, q_init = qInit, maxIt=200, tol = 1e-8)

        qEnd = qSol[0:8]

        if draw:
            self.drawTrajectoryVTK(qInit, qEnd, animate = animate, fold = fold)

        return qEnd, poseError, angleError

    def inverseKinematics(self, endMatrix = None, q_init = None, maxIt=100, tol = 1e-6):
        """
        @brief inverseKinematics: Solve the inverse kinematic problem for the robot.

        Given a wanted end point (inside the boundaries of the robot), the solver will find the parameters
        that move the robot to the position of the wanted end with a tolerance.
        To do so, it uses scipy.optimize.fmin_slsqp minization algorithm using as metric the norm of the
        difference of end point matrices.

        @type endMatrix: numpy.array (4x4)
        @param endMatrix: Wanted end transformation. If None, function returns the original position

        @type q_init: numpy.array (Number of Joints of the robot x 1)
        @param q_init: Initial position of the robot

        @type maxIt: int
        @param maxIt: Maximum number of iteration for the solver

        @type tol: float
        @param tol: tolarance for the solver

        @return q_solution: parameters that place the robot with a Transformation Matrix at
        the end matching endMatrix parameter
        """

        if endMatrix is None:
            return numpy.zeros((self.JointNumber))
        else:

            # Definition of metric function
            def metric_error(params):
                # Calculate using params the Awi
                self.moveRobot(params)
                A_current = self.getRobotEnd()
                error = numpy.sqrt(numpy.sum(((A_current - endMatrix).flatten())**2))
                return error
                return error

            if q_init is None:
                q_init = numpy.zeros((self.JointNumber))

            bounds = list()
            for i in range(self.JointNumber):
                joint = self.JointList[i+1]
                bounds.append((joint.minParam,joint.maxParam))


            q_solved = scipy.optimize.fmin_slsqp( func=metric_error,\
                                              x0=q_init,\
                                              bounds=bounds,\
                                              iprint=0,\
                                              acc=tol,\
                                              iter = maxIt)

            self.moveRobot(q_solved)
            A_current = self.getRobotEnd()

            poseError = numpy.sqrt(numpy.sum(((A_current[0:3,3] - endMatrix[0:3,3]).flatten())**2))
            angleError = numpy.rad2deg(numpy.arccos(numpy.dot(A_current[2,0:3],endMatrix[2,0:3])))

            return q_solved, poseError, angleError

    def drawTrajectoryVTK(self, qInit, qEnd, animate = False, trajectoryColor = [0.0,1.0,1.0], fold = False):


        qInit9 = list(qInit)
        qEnd9 = list(qEnd)

        for i in range(2):
            qInit9.append(0)
            qEnd9.append(0)

        if fold:
            # Fold applicator at beginning before move
            qTrajectories = list()

            # Init to fold state
            qStep1 = list(qInit9)
            qStep1[4] = 0
            qStep1[5] = 0
            qStep1[6] = 15

            # move first stage
            qStep2 = list(qStep1)
            qStep2[0] = qEnd9[0]
            qStep2[1] = qEnd9[1]

            # Unfold liac
            qStep3 = list(qStep2)
            qStep3[4] = 0.89
            qStep3[6] = -30

            # Move 2nd step
            qStep4 = list(qStep3)
            qStep4[2] = qEnd9[2]
            qStep4[3] = qEnd9[3]

            # Fold but Updown
            qStep5 = list(qEnd9)
            qStep5[4] = 0.89

            # Updown movement
            # qEnd

            qTrajectories.append(qInit9)
            qTrajectories.append(qStep1)
            qTrajectories.append(qStep2)
            qTrajectories.append(qStep3)
            qTrajectories.append(qStep4)
            qTrajectories.append(qStep5)
            qTrajectories.append(qEnd9)

            self.TrajectoriesLastSolution = qTrajectories.copy()

            self.drawTrajectoryFoldVTK(qTrajectories, animate = animate, trajectoryColor = trajectoryColor)

        else:
            PyRobotics.RobotArm.drawTrajectoryVTK(self, qInit9, qEnd9, animate = animate, trajectoryColor = trajectoryColor)

    def drawTrajectoryFoldVTK(self, qTrajectories, animate = False, trajectoryColor = [1.0,0.0,1.0]):


        self.trajectory = self.createTrajectory(qTrajectories[0], qTrajectories[1], N =20)

        for i in range(1,len(qTrajectories)-1):
            qi = qTrajectories[i]
            qe = qTrajectories[i+1]
            self.trajectory = numpy.concatenate((self.trajectory,self.createTrajectory(qi, qe, N =20)))


        Npoints, Nparams = self.trajectory.shape
        points = numpy.zeros((Npoints,3))

        for i in range(Npoints):
            q = self.trajectory[i,:]
            self.moveRobot(q)
            points[i,:] = self.JointList[-1].Axis[0:3,3]


        pointsVTK = vtk.vtkPoints()
        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(Npoints)

        for i in range(Npoints):
            pointsVTK.InsertNextPoint(points[i,0],points[i,1],points[i,2])
            polyLine.GetPointIds().SetId(i,i)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyLine)
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(pointsVTK)
        polyData.SetLines(cells)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)
        trajectoryActor = vtk.vtkActor()
        trajectoryActor.SetMapper(mapper)
        trajectoryActor.GetProperty().SetColor(trajectoryColor[0],trajectoryColor[1],trajectoryColor[2])


        # Init to first point
        q = self.trajectory[0,:]
        self.moveRobot(q)

        # Update robot segments and models
        self.JointList[-1].update()

        # Add trajectory to axis
        if animate:
            self.timer_count = 0
            self.timer_countMax = Npoints
            self.timerBool = 1.0

        self.drawVTK(drawAxis = True, actorsToAdd = [trajectoryActor], animate = animate)

    def drawVTK(self, drawAxis = True, actorsToAdd = None, animate = False):

        if actorsToAdd is None:
            actorsToAdd = list()

        if self.calibrated:
            actorsToAdd.append(self.Opticalaxes)

        if self.objetiveFixed:
            actorsToAdd.append(self.ApplicatorObjetiveaxes)
            actorsToAdd.append(self.modelActorApplicator)

        PyRobotics.RobotArm.drawVTK(self, drawAxis = drawAxis, actorsToAdd = actorsToAdd, animate = animate)

    def moveRobot(self, params):

        PyRobotics.RobotArm.moveRobot(self, params = params)

        # Update Optical Axis
        if self.calibrated:

            self.TW2Op = numpy.dot(self.Joint4.Axis, self.TJoint4ToOptical)

            ## Add the Optical Cameras Axis
            self.TOptical.SetMatrix([self.TW2Op[0,0],self.TW2Op[0,1],self.TW2Op[0,2],self.TW2Op[0,3],
                            self.TW2Op[1,0],self.TW2Op[1,1],self.TW2Op[1,2],self.TW2Op[1,3],
                            self.TW2Op[2,0],self.TW2Op[2,1],self.TW2Op[2,2],self.TW2Op[2,3],
                            self.TW2Op[3,0],self.TW2Op[3,1],self.TW2Op[3,2],self.TW2Op[3,3]])
            self.Opticalaxes.SetUserTransform(self.TOptical)




