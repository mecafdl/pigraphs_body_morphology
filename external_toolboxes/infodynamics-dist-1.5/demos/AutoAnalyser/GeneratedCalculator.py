from jpype import *
import numpy
# Our python data file readers are a bit of a hack, python users will do better on this:
sys.path.append("/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/demos/python")
import readFloatsFile

# Add JIDT jar library to the path
jarLocation = "/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# 0. Load/prepare the data:
dataRaw = readFloatsFile.readFloatsFile("/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/demos/data/panda_data.txt")
# As numpy array:
data = numpy.array(dataRaw)
# 1. Construct the calculator:
calcClass = JPackage("infodynamics.measures.continuous.kernel").MutualInfoCalculatorMultiVariateKernel
calc = calcClass()
# 2. Set any properties to non-default values:
calc.setProperty("KERNEL_WIDTH", "0.8")

# Compute for all pairs:
for s in range(42):
    for d in range(42):
        # For each source-dest pair:
        if (s == d):
            continue
        source = data[:, s]
        destination = data[:, d]

        # 3. Initialise the calculator for (re-)use:
        calc.initialise()
        # 4. Supply the sample data:
        calc.setObservations(source, destination)
        # 5. Compute the estimate:
        result = calc.computeAverageLocalOfObservations()

        print("MI_Kernel(col_%d -> col_%d) = %.4f bits" %
            (s, d, result))
