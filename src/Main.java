import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import java.util.Arrays;


import java.util.Random;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int rows = 9;
        int columns = 7;
        Array2DRowRealMatrix document = new Array2DRowRealMatrix(rows, columns);
        Array2DRowRealMatrix reducedMatrix;
        EuclideanDistance euclideanDistance = new EuclideanDistance();
        int counter = 0;
        // Fill document
        int value = 0;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                value = (int) (Math.floor(Math.random() * 16));
                document.setEntry(i, j, value);
                counter++;
            }
            //System.out.println("-------------------------------");
        }
        //System.out.println(counter);
        //System.out.println("--------------------------------------");
        SingularValueDecomposition matrix = new SingularValueDecomposition(document);
        reducedMatrix = (Array2DRowRealMatrix) matrix.getVT();
        System.out.println("\n" + "Reduced Matrix ");
        printResult(reducedMatrix);

        int k = getK(matrix);
        double[][] dStar = eraseRows(reducedMatrix, k, columns);
        System.out.println("dstar lengths");
        System.out.println(dStar.length);
        System.out.println(dStar[0].length);
        System.out.println("-----------");
        int stop = 0;
        while (stop != -1) {
            System.out.println("What documents do you want to compare?");
            int doc1Column = scanner.nextInt();
            int doc2Column = scanner.nextInt();
            double[] doc1 = getColumnValues(dStar, doc1Column);
            double[] doc2 = getColumnValues(dStar, doc2Column);

            double ed = euclideanDistance.compute(doc1, doc2);
            System.out.println("Euclidean distance " + ed);
            System.out.println("Enter -1 to stop comparing, 0 or another number to compare other two documents.");
            stop = scanner.nextInt();
        }

        //SECOND PART

        double[] biggestValues = biggerValues(document, rows, columns);
        double[] q = new double[]{1, 0, 1, 1, 0, 1, 1, 1, 0};
        double[] qPrime = new double[rows];
        System.out.println("biggestValues");
        for (double temp : biggestValues) {
            System.out.println(temp);
        }
        System.out.println("-----------------------------------");
        for (int i = 0; i < rows; i++) {
            qPrime[i] = biggestValues[i] * q[i];
        }
        System.out.println("Lenghts of ts");
        Array2DRowRealMatrix uMatrix = (Array2DRowRealMatrix) matrix.getU();
        double[][] uMatrixStar = eraseColumns(uMatrix, k, rows);
        System.out.println("umatrix = " + uMatrixStar.length + " " + uMatrixStar[0].length);
        Array2DRowRealMatrix sMatrix = (Array2DRowRealMatrix) matrix.getS();
        double[][] sMatrixStar = computeSMatrixStar(sMatrix, k);
        System.out.println("smatrix = " + sMatrixStar.length + " " + sMatrixStar[0].length);
        System.out.println("-------------------------------------------------");


        double[][] uSMatrixStar = multiplyArray2DMatrices(uMatrixStar, sMatrixStar);
        /*System.out.println(uSMatrixStar.length);
        System.out.println(uSMatrixStar[0].length);*/


        System.out.println("Q Prime");
        for (double temp : qPrime) {
            System.out.println(temp);
        }


        System.out.println(uSMatrixStar.length); // 7
        System.out.println(uSMatrixStar[0].length); // 2

        double[][] usMatrixTransposed = transpose(uSMatrixStar);
        System.out.println(usMatrixTransposed.length); // 2
        System.out.println(usMatrixTransposed[0].length); // 7

        double[][] qStar = new double[usMatrixTransposed.length][1];
        //for (int i = 0; i < qPrime.length; i++) {
        for (int j = 0; j < usMatrixTransposed.length; j++) {
            qStar[j][0] = 0;
            for (int l = 0; l < qPrime.length; l++) {
                qStar[j][0] += qPrime[l] * usMatrixTransposed[j][l];
            }
        }
        //}
        double[] cleanedQStar = new double[usMatrixTransposed.length];
        System.out.println("Qstar values");
        for (int i = 0; i < usMatrixTransposed.length; i++) {
            System.out.println(qStar[i][0]);
            cleanedQStar[i] = qStar[i][0];
        }


       /* System.out.println("Qstar: " + "\n");
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < columns; j++) {
                cleanedQStar[j] = qStar[i][j];
                System.out.println(qStar[i][j]);
            }
        }*/


        System.out.println("Doc comparison cleanedQstar");
        //Compare qStar with d*T (reduced matrix) docs
        double[] temporalDocument = new double[columns];
        double temporalDistance = 0;
        double[] distancesArray = new double[columns];
        double[] sortedDistancesArray = new double[columns];
        for (int i = 0; i < columns; i++) {
            temporalDocument = getColumnValues(dStar, i);
            temporalDistance = euclideanDistance.compute(temporalDocument, cleanedQStar);
            distancesArray[i] = temporalDistance;
            sortedDistancesArray[i] = temporalDistance;
        }
        Arrays.sort(sortedDistancesArray);
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < columns; j++) {
                if (distancesArray[i] == sortedDistancesArray[j]) {
                    System.out.println("Document number " + (i + 1) + " is the " + (j + 1) + " closest one, and the distance is " + sortedDistancesArray[j]);
                }
            }

        }
        System.out.println("----------------------------------------");


        //System.out.println(getK(matrix));

        //printResult(reducedMatrix);
    }


    public static void printResult(Array2DRowRealMatrix matrix) {
        double[][] result = matrix.getData();
        int counter = 0;
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                counter++;
                // Print result
                //System.out.println(result[i][j]);
            }
        }
        //System.out.println(counter);

    }


    public static double[] getColumnValues(double[][] matrix, int column) {
        //double[][] data = matrix.getData();
        double[] result = new double[matrix.length];

        //System.out.println("Values found in each column");

        for (int i = 0; i < matrix.length; i++) {
            result[i] = matrix[i][column];
            // System.out.println(result[i]);
        }
        //System.out.println("----------------------------------------------");
        return result;
    }


    public static double[][] eraseRows(Array2DRowRealMatrix matrix, int k, int columns) {
        double[][] data = matrix.getData();
        double[][] result = new double[k + 1][columns];

        for (int i = 0; i <= k; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = data[i][j];
            }
        }
        return result;
    }


    public static double[][] eraseColumns(Array2DRowRealMatrix matrix, int k, int rows) {
        double[][] data = matrix.getData();
        double[][] result = new double[rows][k + 1];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < k + 1; j++) {
                result[i][j] = data[i][j];
            }
        }
        return result;
    }


    public static double[][] computeSMatrixStar(Array2DRowRealMatrix matrix, int k) {
        double[][] data = matrix.getData();
        double[][] result = new double[k + 1][k + 1];

        for (int i = 0; i < k + 1; i++) {
            for (int j = 0; j < k + 1; j++) {
                result[i][j] = data[i][j];
            }
        }
        return result;
    }


    public static int getK(SingularValueDecomposition matrix) {
        double[] singularValues = matrix.getSingularValues();
        double limit = singularValues[1] * 0.05;
        int k = 1;

        for (int i = 1; i < singularValues.length - 1; i++) {
            if (singularValues[i] - singularValues[i + 1] > limit) {
                return k;
            } else {
                k++;
            }
        }
        return k;
    }


    public static double[] biggerValues(Array2DRowRealMatrix document, int rows, int columns) {
        double[][] matrix = document.getData();
        double[] result = new double[rows];
        double maxRowValue = 0;

        //System.out.println("Bigger values");
        for (int i = 0; i < rows; i++) {
            //  System.out.println("i=" + i);
            for (int j = 0; j < columns; j++) {
                //    System.out.println("Value = " + matrix[i][j]);
                if (matrix[i][j] > maxRowValue) {
                    //      System.out.println("biggest value = "+ matrix[i][j]);
                    result[i] = matrix[i][j];
                    maxRowValue = matrix[i][j];
                }
            }
            maxRowValue = 0;
            //System.out.println("-----------------------------------");
        }

        return result;
    }


    public static double[][] multiplyArray2DMatrices(double[][] firstMatrix, double[][] secondMatrix) {
        System.out.println("MULTIPLICATION");
        System.out.println("Matrix1 rxc = " + firstMatrix.length + " " + firstMatrix[0].length);
        System.out.println("Matrix2 rxc = " + secondMatrix.length + " " + secondMatrix[0].length);
        double[][] result = new double[firstMatrix.length][secondMatrix[0].length];

        for (int i = 0; i < firstMatrix.length; i++) {
            for (int j = 0; j < secondMatrix[0].length; j++) {
                result[i][j] = 0;
                for (int k = 0; k < firstMatrix[0].length; k++) {
                    result[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                }
            }
        }
        System.out.println("--------------------------------------------");
        return result;
    }


    public static double[][] transpose(double[][] matrix) {
        double[][] result = new double[matrix[0].length][matrix.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                result[j][i] = matrix[i][j];
            }
        }

        return result;
    }

}
