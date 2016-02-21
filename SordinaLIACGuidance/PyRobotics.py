# import plotly
import numpy
import scipy
import copy
import os
import vtk

def DHtoT(A, alpha, D, offset):
    """
    @brief DHtoT Return the transformation matrix (4x4) from the Denavit-Hartenberg parameters

    @type A: double
    @param A: Axis X displacement (fixed parameter for the joint definition)

    @type alpha: double
    @param alpha: Axis X rotation (degrees) (fixed parameter for the joint definition)

    @type D: double
    @param D: Axis Z displacement (this parameter could be variable if the joint is periscopic)

    @type offset: double
    @param offset: Axis Z rotation (degrees) (this parameter could be variable if the joint is rotational)

    @rtype: numpy.array
    @return: Returns a numpy.array (4x4) which is the transformation from the Denavit-Hartenberg parameters
    """
    # A, alpha, D, offset.
    # Tx, Rx, Tz, Rz.

    # A = Rz * Tz * Tx * Rx


    # X Translation
    Tx = numpy.eye(4)
    Tx[0,3] = A

    # X Rotation
    Rx = numpy.eye(4)
    cx = numpy.cos(numpy.deg2rad(alpha))
    sx = numpy.sin(numpy.deg2rad(alpha))
    Rx[1,1] = cx
    Rx[1,2] = -sx
    Rx[2,1] = sx
    Rx[2,2] = cx

    # Z Rotation
    Rz = numpy.eye(4)
    cz = numpy.cos(numpy.deg2rad(offset))
    sz = numpy.sin(numpy.deg2rad(offset))
    Rz[0,0] = cz
    Rz[0,1] = -sz
    Rz[1,0] = sz
    Rz[1,1] = cz

    # Z Translation
    Tz = numpy.eye(4)
    Tz[2,3] = D

    # Premultiplication
    T = numpy.dot(Rz,numpy.dot(Tz,numpy.dot(Tx,Rx)))
    return T

def createAxisFigureObject(T,label = '', s = 1.0):
    """
    @brief createAxisFigureObject Create object of plotly for draw an axis according to transformation T (4x4)

    @type T: numpy.array
    @param T: Transformation that defines the rotation and center of the axis

    @return: List of drawable objects
    """

    # Create point for center
    center = plotly.graph_objs.Scatter3d(
        x = numpy.dot(T,numpy.array([0,0,0,1]))[0],
        y = numpy.dot(T,numpy.array([0,0,0,1]))[1],
        z = numpy.dot(T,numpy.array([0,0,0,1]))[2],
        mode='markers',
        marker=dict(
            color='black',
            size=12,
            symbol='circle',
            line=dict(
                color='rgb(204, 204, 204)',
                width=1
            ),
            opacity=0.9
        )
    )

    # create Axis X
    axisX = plotly.graph_objs.Scatter3d(
                name = "X " + label,
                x = [numpy.dot(T,numpy.array([0,0,0,1]))[0], numpy.dot(T,numpy.array([s,0,0,1]))[0]],
                y = [numpy.dot(T,numpy.array([0,0,0,1]))[1], numpy.dot(T,numpy.array([s,0,0,1]))[1]],
                z = [numpy.dot(T,numpy.array([0,0,0,1]))[2], numpy.dot(T,numpy.array([s,0,0,1]))[2]],
                mode='lines+text',
                line=plotly.graph_objs.Line(
                    color='red',
                    width=4
                ),
                text=['', 'X' + label],
                textposition='top',
                textfont=dict(
                    size=10,
                    color='red'
                )
            )

    # Create Axis Y
    axisY = plotly.graph_objs.Scatter3d(
                name = "Y " + label,
                x = [numpy.dot(T,numpy.array([0,0,0,1]))[0], numpy.dot(T,numpy.array([0,s,0,1]))[0]],
                y = [numpy.dot(T,numpy.array([0,0,0,1]))[1], numpy.dot(T,numpy.array([0,s,0,1]))[1]],
                z = [numpy.dot(T,numpy.array([0,0,0,1]))[2], numpy.dot(T,numpy.array([0,s,0,1]))[2]],
                mode='lines+text',
                line=plotly.graph_objs.Line(
                    color='green',
                    width=4
                ),
                text=['', 'Y' + label],
                textposition='top',
                textfont=dict(
                    size=10,
                    color='green'
                )
            )
    # Create Axis Z
    axisZ = plotly.graph_objs.Scatter3d(
                name = "Z " + label,
                x = [numpy.dot(T,numpy.array([0,0,0,1]))[0], numpy.dot(T,numpy.array([0,0,s,1]))[0]],
                y = [numpy.dot(T,numpy.array([0,0,0,1]))[1], numpy.dot(T,numpy.array([0,0,s,1]))[1]],
                z = [numpy.dot(T,numpy.array([0,0,0,1]))[2], numpy.dot(T,numpy.array([0,0,s,1]))[2]],
                mode='lines+text',
                line=plotly.graph_objs.Line(
                    color='blue',
                    width=4
                ),
                text=['', 'Z' + label],
                textposition='left',
                textfont=dict(
                    size=10,
                    color='blue'
                )
            )

    data = [center, axisX, axisY, axisZ]
    return data

def draw3DObjectsPyplot(data, title= ''):
    """
    @brief draw3DObjectsPyplot Draw object in list data using plotly

    @type data: list
    @param data: drawable list of plotly object for drawing

    @return: None
    """
    layout = dict(
                title = title,
                autosize=True,
                showlegend=False,
                scene = dict(
                    aspectmode = 'data'
                )
            )

    fig = plotly.graph_objs.Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)

def createAxisZ(T,alphaV = 0.5,length = 1.0):
    """
    @brief createAxisZ Create object of plotly for draw an axis according to transformation T (4x4)

    @type T: numpy.array
    @param T: Transformation that defines the rotation and center of the axis

    @return: List of drawable objects
    """

    f = length/2.0

    # Create Axis Z
    axisZ = plotly.graph_objs.Scatter3d(
                name = "line",
                x = [numpy.dot(T,numpy.array([0,0,-f,1]))[0], numpy.dot(T,numpy.array([0,0,f,1]))[0]],
                y = [numpy.dot(T,numpy.array([0,0,-f,1]))[1], numpy.dot(T,numpy.array([0,0,f,1]))[1]],
                z = [numpy.dot(T,numpy.array([0,0,-f,1]))[2], numpy.dot(T,numpy.array([0,0,f,1]))[2]],
                mode='lines',
                line=plotly.graph_objs.Line(
                    color=u'rgba(0,0,255,' + str(alphaV) + ')',
                    width=1,
                )
            );
    return axisZ

class Joint:
    """@brief Class representing a Joint of a robot

    Initialization is made from Denavit-Hartenberg parameters of the joint inside a chain
    """

    def __init__(self, A, alpha, D, offset, freedom, idx = 0):
        """
        @brief __init__ Create Joint object from Denavit-Hartenberg parameters

        @type A: float
        @param A: Translation in X axis

        @type alpha: float
        @param alpha: Rotation in X axis (degrees)

        @type D: float
        @param D: Translation in Z axis

        @type offset: float
        @param offset: Rotation in Z axis (degrees)

        @type freedom: string
        @param freedom: 'rotation' if is a rotation joint, or 'periscopic' if the joint is periscopic

        @type idx: int
        @param idx: Joint Identificator
        """

        self.A = A
        self.alpha = alpha
        self.D = D
        self.offset = offset
        self.idx = idx

        if freedom == 'rotation':
            self.jointType = 'rotation'
        elif freedom == 'periscopic':
            self.jointType = 'periscopic'

        # degree of freedom of the joint
        self.param = 0.0
        # boundaries of the param
        self.minParam = -numpy.inf
        self.maxParam = numpy.inf
        # True if the joint is the first in a chain
        self.Isfirst = True
        # Transformation matrix from previous joint to current joint
        self.Ai = DHtoT(self.A, self.alpha, self.D, self.offset)
        # Transformation matrix from world coordinates to joint frame
        self.Awi = self.Ai.copy()
        self.Axis = self.Ai.copy()

        # VTK axes
        self.axesActor = vtk.vtkAxesActor()
        self.vtkT = vtk.vtkTransform()
        self.axesActor.SetXAxisLabelText('X ' + str(self.idx))
        self.axesActor.SetYAxisLabelText('Y ' + str(self.idx))
        self.axesActor.SetZAxisLabelText('Z ' + str(self.idx))

        # VTK objects Model
        self.pathToSTLModel = None
        self.STLReader = None
        self.mapper = None
        self.modelActor = None

        self.setModel()


    def setFrame(self,jointFrame = None):
        """
        @brief setFrame Set the previous frame to which current frame is connected

        @type jointFrame: Joint
        @param jointFrame: Joint object defining the previous joint
        """
        if jointFrame is None:
            # The Joint is supposed to be connected to World center
            self.frame = Joint(0,0,0,0,'rotation')
            self.Isfirst = True
        else:
            # The Joint is connected to other frame
            self.frame = jointFrame
            self.Isfirst = False
            #print "Connected {0:d} -> {1:d}".format(self.frame.idx, self.idx)

        self.update()

    def update(self):
        """
        @brief update Update the matrices self.Ai and self.Awi according to params in
        current joint and previous connected joints in the chain
        """

        # This function update new location of axis
        self.Ai = DHtoT(self.A, self.alpha, self.D, self.offset)

        # Calculate new Ai using parameter
        if self.jointType == 'periscopic':
            self.Ai = numpy.dot(self.Ai,DHtoT(0, 0, self.param, 0))
        else:
            self.Ai = numpy.dot(self.Ai,DHtoT(0,0,0, self.param))

        # Calculate new axis
        # Update connected if is part of a chain
        if self.Isfirst:
            self.Awi = self.Ai.copy()
        else:
            self.frame.update()
            self.Awi = numpy.dot(self.frame.Awi, self.Ai)

        self.Axis = self.Awi.copy()

        self.vtkT.SetMatrix([self.Axis[0,0], self.Axis[0,1], self.Axis[0,2], self.Axis[0,3],
                             self.Axis[1,0], self.Axis[1,1], self.Axis[1,2], self.Axis[1,3],
                             self.Axis[2,0], self.Axis[2,1], self.Axis[2,2], self.Axis[2,3],
                             self.Axis[3,0], self.Axis[3,1], self.Axis[3,2], self.Axis[3,3]])

        self.axesActor.SetUserTransform(self.vtkT)

        if self.modelActor is not None:
            self.modelActor.SetUserTransform(self.vtkT)

    def rotate(self, degrees):
        """
        @brief rotate Rotate joint if is a rotation joint.

        @type degrees: float
        @param degrees: angle to rotate from original position of the joint (degrees)
        """
        if self.jointType == 'rotation':
            if (degrees <= self.maxParam) and (degrees >= self.minParam):
                self.param = degrees
                self.update()
            else:
                print "Movement is out of limits"
        elif self.jointType == 'periscopic':
            print "No rotation applied, this joint is not a rotation joint"

    def translate(self, periscopic):
        """
        @brief translate Translate joint if is a periscopic joint.

        @type periscopic: float
        @param periscopic: distance to move the joint from original position
        """
        if self.jointType == 'periscopic':
            if (periscopic <= self.maxParam) and (periscopic >= self.minParam):
                self.param = distance
                self.update()
            else:
                print "Movement is out of limits"
        elif self.jointType == 'rotation':
            print "No translation applied, this joint is not a periscopic joint"

    def move(self, q):
        """
        @brief move: Translate or rotate joint if is a periscopic or rotation joint.

        @type q: float
        @param q: distance/angle to move the joint from original position m/degrees
        """
        self.param = q
        self.update()

    def setLimits(self, minParam, maxParam):
        """
        @brief setLimits: set the limits of the joints

        @type minParam: float
        @param minParam: minimum distance/angle to move the joint from original position m/degrees

        @type maxParam: float
        @param maxParam: maximum distance/angle to move the joint from original position m/degrees
        """
        self.minParam = minParam
        self.maxParam = maxParam

    def setModel(self, pathToSTLModel = '', color = [0,1,1], opacity = 0.4, length = 1.0, radius = 0.1):

        if os.path.isfile(pathToSTLModel):

            self.pathToSTLModel = pathToSTLModel

            self.STLReader = vtk.vtkSTLReader()
            self.STLReader.SetFileName(self.pathToSTLModel)
            self.STLReader.Update()

            self.mapper = vtk.vtkPolyDataMapper()
            self.mapper.SetInputConnection(self.STLReader.GetOutputPort())

            self.modelActor = vtk.vtkActor()
            self.modelActor.SetMapper(self.mapper)

            self.modelActor.GetProperty().SetColor(color[0], color[1], color[2])
            self.modelActor.GetProperty().SetOpacity(opacity)

        else:

            self.tube = vtk.vtkLineSource()
            self.tube.SetPoint1( 0, 0, 0 )
            self.tube.SetPoint2( 0, 0, length )
            self.tube.SetResolution( 100 )

            self.vtkTubeFilter = vtk.vtkTubeFilter()
            self.vtkTubeFilter.SetInputConnection(self.tube.GetOutputPort())
            self.vtkTubeFilter.SetRadius(radius)
            self.vtkTubeFilter.SetNumberOfSides(20)
            self.vtkTubeFilter.CappingOn()

            self.mapper = vtk.vtkPolyDataMapper()
            self.mapper.SetInputConnection(self.vtkTubeFilter.GetOutputPort())

            self.modelActor = vtk.vtkActor()
            self.modelActor.SetMapper(self.mapper)

            self.modelActor.GetProperty().SetColor(color[0], color[1], color[2])
            self.modelActor.GetProperty().SetOpacity(opacity)

    def setAxesLength(self, axesLength = 0.2):

        self.axesActor.SetTotalLength(axesLength,axesLength,axesLength)

    def __repr__(self):
        self.update()
        s = "Joint number [{0:s}]".format(str(self.idx))
        s = s + '\n' + "    - Type: " + self.jointType
        s = s + '\n' + "    - Parameters (A, alpha, D, offset) = ({0:+.2f}, {1:+.2f}, {2:+.2f}, {3:+.2f})".format(self.A,\
                                                                                                         self.alpha,\
                                                                                                         self.D,\
                                                                                                         self.offset)
        s = s + '\n' + "    - Center of joint [{0:+.2f}, {1:+.2f}, {2:+.2f}]".format(self.Awi[0,3],\
                                                                            self.Awi[1,3],\
                                                                            self.Awi[2,3])

        if self.pathToSTLModel is not None:
            s = s + '\n' + "    - Associated Model: " + self.pathToSTLModel
        else:
            s = s + '\n' + "    - Associated Model: VTK cylinder"
        return s

    def draw(self):

        # Update location and movement
        self.update()

        # create a rendering window and renderer
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.SetSize(1200, 900)
        renWin.AddRenderer(ren)

        # Create camera
        camera = vtk.vtkCamera()
        #camera.SetPosition(0.5,10.0,0.0);
        #camera.Roll(-90)
        #ren.SetActiveCamera(camera)

        # create a renderwindowinteractor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # Add actors
        ren.AddActor(self.modelActor)
        ren.AddActor(self.axesActor)

        # enable user interface interactor
        iren.Initialize()
        renWin.Render()
        style = vtk.vtkInteractorStyleMultiTouchCamera()
        iren.SetInteractorStyle(style)
        iren.Start()

class RobotArm:
    """@brief Class representing a chain of joints of a robot arm

    Initialization is made from Denavit-Hartenberg parameters of the joint inside a chain
    """

    def __init__(self):
        """
            @brief create a robotArm
        """
        # List of joints of the chaing
        self.JointList = list()
        # Create the first joint as World
        self.JointList.append(Joint(0,   0, 0,   0, 'rotation', idx = 'World'))
        self.JointNumber = 0
        self.MaximumRange = None

    def moveRobot(self, params):
        """
        @brief moveRobot: move the robot using the vector params

        @type params: numpy.array 1D
        @param params: parameters of the joints to be moved. The order start in the first joint that is connected to the base

        """
        params = numpy.array(params)
        paramsX = numpy.concatenate((numpy.zeros((1)),params))
        N = len(paramsX)
        if N == len(self.JointList):
            for (j,q) in zip(self.JointList,paramsX):
                j.move(q)
        else:
            print "Please, check the number of parameters passed to the moveRobot() function"

    def connectJoint(self, newJoint):
        """
        @brief connectJoint: add a Joint to the end of the chain

        @type newJoint: Joint
        @param newJoint: new Joint to be added, already initialized
        """
        self.JointNumber += 1
        if len(self.JointList) > 0:
            newJoint.setFrame(self.JointList[-1])

        # Reset Joint
        newJoint.move(0.0)
        self.JointList.append(newJoint)

    def resetRobot(self):
        """
        @brief resetRobot: clear chain and delete joints
        """
        self.JointList = list()
        self.JointList(Joint(0,   0, 0,   0, 'rotation', idx = 'World'))

    def getRobotEnd(self):
        """
        @brief getRobotEnd: return the transformation matrix that moves the base origin to the end of the robot chain

        @return: numpy.array((4,4)) matrix transformation

        """

        return self.JointList[-1].Awi

    def getJoint(self, n):
        """
        @brief getJoint: returns the joint at location n into the chain

        @type n: int
        @param n: 0-index of the joint into the robot chain

        @return Joint at this location into the chain order
        """
        return self.JointList[n+1]

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
            return numpy.zeros((self.JointNumber)), 0, 0
        else:

            # Definition of metric function
            def metric_error(params):
                # Calculate using params the Awi
                self.moveRobot(params)
                A_current = self.getRobotEnd()
                error = numpy.sqrt(numpy.sum(((A_current - endMatrix).flatten())**2))
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

            error = numpy.sqrt(numpy.sum(((A_current - endMatrix).flatten())**2))

            return q_solved, error

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

            for i in range(self.JointNumber - 1): # Last one is not rotable (fake)
                joint = self.getJoint(i)
                print "Joint", joint.idx, "Max", joint.maxParam, "Min", joint.minParam
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
        return self.MaximumRange

    def draw(self, title = '', s =1.0, drawRange = False, drawAll=True, data = None, drawAxis = True):
        """
        @brief draw: crate the plotly offline figure for the robot

        @type title: string
        @param title: Title of the figure

        @type drawRange: Bool
        @param drawRange: if true, the maximumRange will be computed and added to the figure

        @type drawAll: Bool
        @param drawAll: if true, the figure is created, else the figure-object is only returned with no drawing

        @return drawable list of the required objects in plotly configuration
        """

        if data is None:
            data = list()
        AxisObjects = list()

        centers = numpy.zeros((len(self.JointList),3))
        for i, j in enumerate(self.JointList):
            AxisObjects.append(createAxisFigureObject(j.Axis, str(j.idx), s = s))
            centers[i,:] = numpy.dot(j.Axis,numpy.array([0,0,0,1]))[0:3]


        robotLine = plotly.graph_objs.Scatter3d(
                    name = "RobotLine ",
                    x = centers[:,0],
                    y = centers[:,1],
                    z = centers[:,2],
                    mode='lines',
                    line=plotly.graph_objs.Line(
                        color='black',
                        width=8
                    )
                )


        if drawAxis:
            for ax in AxisObjects:
                for element in ax:
                    data.append(element)

        data.append(robotLine)
        if drawRange:
            data.append(self.computeMaximumRange())
        if drawAll:
            draw3DObjectsPyplot(data, title = title)
        return data

    def drawTrajectory(self, qmatrix):
        Npoints, Nparams = qmatrix.shape
        points = numpy.zeros((Npoints,3))

        for i in range(Npoints):
            q = qmatrix[i,:]
            self.moveRobot(q)
            points[i,:]= self.getRobotEnd()[0:3,3]

            # Create Axis Z
        pointsFigure = plotly.graph_objs.Scatter3d(
                            name = "line",
                            x = points[:,0],
                            y = points[:,1],
                            z = points[:,2],
                            mode='markers',
                            marker=dict(
                                        size=2,
                                        line=dict(
                                            color='rgba(217, 217, 217, 0.14)',
                                            width=0.5
                                        ),
                                        opacity=0.8
                                    )
                        )
        #self.draw(title = 'Trajectory', drawRange = False, drawAll=True, data = [pointsFigure])
        return pointsFigure, points

    def drawVTK(self, drawAxis = True, actorsToAdd = None, animate = False):

        # Update joints
        self.JointList[-1].update()

        # create a rendering window and renderer
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.SetSize(1200, 900)
        renWin.AddRenderer(ren)

        # Create camera
        camera = vtk.vtkCamera()
        #camera.SetPosition(0.5,10.0,0.0);
        #camera.Roll(-90)
        #ren.SetActiveCamera(camera)

        # create a renderwindowinteractor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # Add actors

        for j in self.JointList[1::]:
            ren.AddActor(j.modelActor)
            if drawAxis:
                ren.AddActor(j.axesActor)

        # Axes for world
        ren.AddActor(self.JointList[0].axesActor)

        # Added Actors
        if actorsToAdd is not None:
            for a in actorsToAdd:
                ren.AddActor(a)

        # enable user interface interactor
        iren.Initialize()
        renWin.Render()

        if animate:
            iren.AddObserver('TimerEvent', self.animationCallBack)
            timerId = iren.CreateRepeatingTimer(100);


        style = vtk.vtkInteractorStyleMultiTouchCamera()
        iren.SetInteractorStyle(style)
        iren.Start()

    def createTrajectory(self,qInit,qEnd, N = 10):


        NumberOfMovements = 0
        InvalidMovements = list()
        for i in range(self.JointNumber):
            if qInit[i] != qEnd[i]:
                NumberOfMovements = NumberOfMovements + 1
            else:
                InvalidMovements.append(i)

        qmatrix = numpy.zeros((N * NumberOfMovements,self.JointNumber))

        for i in range(self.JointNumber):
            qmatrix[:,i] = qInit[i]

        j = 0
        for i in range(self.JointNumber):
            if qInit[i] != qEnd[i]:
                qmatrix[j * N:(N + j * N), i] = numpy.linspace(qInit[i],qEnd[i],N)
                qmatrix[(N + j * N)::, i] = qEnd[i]
                j = j + 1

        return qmatrix

    def drawTrajectoryVTK(self, qInit, qEnd, animate = False, trajectoryColor = [1.0,0.0,1.0]):

        self.trajectory = self.createTrajectory(qInit, qEnd, N =20)

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

    def animationCallBack(self, obj, event):

        # Update timer count
        self.timer_count = self.timer_count % self.timer_countMax
        if self.timer_count == self.timer_countMax - 1:
            self.timerBool = self.timerBool * (-1)

        # Get the trajectory point
        if self.timerBool > 0:
            q = self.trajectory[self.timer_count,:]

        else:
            q = self.trajectory[self.timer_countMax-1,:]
            if self.timer_count == numpy.round(self.timer_countMax*0.2):
                self.timer_count = -1
                self.timerBool = 1

        self.moveRobot(q)
        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1








