package infodynamics.demos.autoanalysis;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

import infodynamics.measures.continuous.*;
import infodynamics.measures.continuous.kernel.*;

public class GeneratedCalculator {

  public static void main(String[] args) throws Exception {

    // 0. Load/prepare the data:
    String dataFile = "/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/demos/data/panda_data.txt";
    ArrayFileReader afr = new ArrayFileReader(dataFile);
    double[][] data = afr.getDouble2DMatrix();
    // 1. Construct the calculator:
    MutualInfoCalculatorMultiVariateKernel calc;
    calc = new MutualInfoCalculatorMultiVariateKernel();
    // 2. Set any properties to non-default values:
    calc.setProperty(MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME,
        "0.8");
    
    // Compute for all pairs:
    for (int s = 0; s < 42; s++) {
        for (int d = 0; d < 42; d++) {
            // For each source-dest pair:
            if (s == d) {
                continue;
            }
            double[] source = MatrixUtils.selectColumn(data, s);
            double[] destination = MatrixUtils.selectColumn(data, d);

            // 3. Initialise the calculator for (re-)use:
            calc.initialise();
            // 4. Supply the sample data:
            calc.setObservations(source, destination);
            // 5. Compute the estimate:
            double result = calc.computeAverageLocalOfObservations();

            System.out.printf("MI_Kernel(col_%d -> col_%d) = %.4f bits\n",
                s, d, result);
        }
    }
  }
}

