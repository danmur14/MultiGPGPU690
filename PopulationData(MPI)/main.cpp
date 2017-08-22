#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <mpi.h>

std::string MPIstrategy, MPIoperation, MPIrelation; // (sr or bg), (max, min, avg, number), (gt, lt) the args passed when running the program
std::string columnLetters[] =
    {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
     "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ",
     "BA", "BB", "BC", "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BK", "BL", "BM", "BN", "BO", "BP", "BQ", "BR", "BS", "BT", "BU", "BV", "BW", "BX", "BY", "BZ",
     "CA", "CB", "CC", "CD", "CE", "CF", "CG", "CH", "CI", "CJ", "CK", "CL", "CM", "CN", "CO", "CP", "BQ", "CR", "CS", "CT", "CU", "CV", "CW", "CX", "CY", "CZ",
     "DA", "DB", "DC", "DD", "DE", "DF", "DG", "DH", "DI", "DJ", "DK", "DL", "DM"}; // an easy way to map the letters on the columns to column indexes

int main(int argc, char **argv)
{
    // MPI calls for world
    MPI_Init(NULL, NULL);
    // returns what rank this process is
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // number of processes in MPI world
    int n_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes); //check how many processes we have

    MPIstrategy = argv[1];  // strategy to use from command arguments, sr or bg
    MPIoperation = argv[2]; //operation used on the data, max min avg number

    // if n_processes doesn't divide 500 evenly, abort
    if (MPIstrategy == "sr" && 500 % n_processes != 0)
    {
        if (rank == 0) //only rank 0 should say error
        {
            std::cout << n_processes << " does not divide 500 evenly. Exiting..." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::map<std::string, int> letterColumnMap; // c++ maps column letters to column indexes using columnLetters array
    std::string columnHeaders[117];             // column header names
    std::string states[500][2];                 // array for name and abbreviations of states

    std::string row;        // string for each row of the csv file when parsing
    double values[500][115]; // data from the csv is placed in this array

    // rank 0 processes will handle parsing of csv file
    if (rank == 0)
    {
        // maps columnLetter to a column index
        for (int i = 0; i < 117; i++)
        {
            letterColumnMap[columnLetters[i]] = i;
        }

        // open the csv file and start parsing the data
        std::ifstream csv;
        csv.open("500_Cities__City-level_Data__GIS_Friendly_Format_.csv");
        if (csv)
        {
            getline(csv, row);
            std::stringstream columnStream(row); // stream for the header row containing column names
            std::string columnName;

            // Read the first line, which contains the columns names
            int columnIndex = 0;
            while (getline(columnStream, columnName, ','))
            {
                columnHeaders[columnIndex] = columnName; // store name in respective column
                columnIndex++;
            }

            // iterate through each row and place data in value array
            int i = 0;                // row index
            while (getline(csv, row)) // not at the end of file
            {
                std::stringstream rowStream(row); // stream for each row containing data
                std::string value;                // value in location

                int j = 0;                             // column index
                while (getline(rowStream, value, ',')) // get value per column on each row
                {
                    bool trueValue = true; // not a ranged value

                    // check if value has quotes and place a 0 instead
                    if (value.find('\"') != std::string::npos)
                    {
                        trueValue = false;
                        values[i][j - 2] = 0; // replace with 0 where ranged and unused
                        j++;
                        getline(rowStream, value, ','); // ranged values have commas, so get the second half before getting the next column value
                    }

                    if (trueValue) //value should be input to the data array
                    {
                        if (j < 2) // first two columns are states and abbreviations
                        {
                            states[i][j] = value;
                        }
                        else // normal data
                        {
                            values[i][j - 2] = std::stof(value);
                        }
                        j++;
                    }
                }
                i++;
            }
            csv.close();
        }
        else
        {
            std::cout << "File could not be opened. Exiting...";
        }
    }

    // MPI_MAXLOC and MPI_MINLOC requires a structure, so send this struct back so rank 0 does not have to calculate the index of desired value, and we do not have to iterate over the data unneccessarily
    struct messageData
    {
        double value; // the value to return from the process, (max, min...)
        int index;   //the index in the total data array
    } in, out;

    // use scatter and reduce
    if (MPIstrategy == "sr")
    {
        int sendCount = 500 / n_processes; //messages per process
        double sendColumn[500], recvBuffer[sendCount];//array to send, results in total and local,

        // rank 0 fills sendColumn with values from appropriate column
        if (rank == 0)
        {
            for (int i = 0; i < 500; i++)
            {
                sendColumn[i] = values[i][letterColumnMap[argv[3]] - 2];
            }
        }

        // rank 0 scatters the filled sendColumn and is received in the recvBuffer
        MPI_Scatter(&sendColumn, sendCount, MPI_DOUBLE, recvBuffer, sendCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // get operation (min, max, avg, number)
        MPIoperation = argv[2];

        if (MPIoperation == "max")
        {
            in.value = recvBuffer[0];
            in.index = 0;

            // iterate on received values to find max
            for (int i = 1; i < sendCount; i++)
            {
                if (recvBuffer[i] > in.value)
                {
                    in.value = recvBuffer[i]; // max value from recvBuffer
                    in.index = i;             // the index of the max value
                }
            }
            in.index = rank * sendCount + in.index; // rank * sendCount is the section received by a particular process, offset by the index

            // send back the max value and the index that it is in
            MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                double totalMax = out.value;
                int maxIndex = out.index;
                std::cout << std::setprecision(15) << states[maxIndex][1] << ", " << states[maxIndex][2] << ", " << columnHeaders[letterColumnMap[argv[3]]] << " = " << totalMax << std::endl;
            }
        }
        else if (MPIoperation == "avg")
        {
            double localSum = 0;
            for (int i = 0; i < sendCount; i++)
            {
                localSum += recvBuffer[i];
            }
            double totalSum;

            MPI_Reduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                std::cout << "Average " << columnHeaders[letterColumnMap[argv[3]]] << " = " << totalSum / (n_processes * sendCount) << std::endl;
            }
        }
        else if (MPIoperation == "min")
        {
            in.value = recvBuffer[0];
            in.index = 0;

            for (int i = 1; i < sendCount; i++)
            {
                if (recvBuffer[i] < in.value)
                {
                    in.value = recvBuffer[i];
                    in.index = i;
                }
            }
            in.index = rank * sendCount + in.index;

            MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                double overallMin = out.value;
                int minIndex = out.index;

                std::cout << states[minIndex][1] << ", " << states[minIndex][2] << ", " << columnHeaders[letterColumnMap[argv[3]]] << " = " << overallMin << std::endl;
            }
        }
        else if (MPIoperation == "number")
        {
            std::string MPIrelation = argv[4];
            double targetValue = std::stof(argv[5]); // convert string to double

            if (MPIrelation == "gt")
            {
                int localGreater = 0;
                for (int i = 0; i < sendCount; i++)
                {
                    if (recvBuffer[i] > targetValue)
                    {
                        localGreater++;
                    }
                }
                int totalGreater;

                MPI_Reduce(&localGreater, &totalGreater, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if (rank == 0)
                {
                    std::cout << "Number cities with " << columnHeaders[letterColumnMap[argv[3]]] << " gt " << targetValue << " = " << totalGreater << std::endl;
                }
            }
            else if (MPIrelation == "lt")
            {
                int localLess = 0;
                for (int i = 0; i < sendCount; i++)
                {
                    if (recvBuffer[i] < targetValue)
                    {
                        localLess++;
                    }
                }
                int totalLess;

                MPI_Reduce(&localLess, &totalLess, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if (rank == 0)
                {
                    std::cout << "Number cities with " << columnHeaders[letterColumnMap[argv[3]]] << " lt " << targetValue << " = " << totalLess << std::endl;
                }
            }
        }
    }
    else if (MPIstrategy == "bg")
    {
        double sendColumns[n_processes][500];

        // rank 0 populates the array with a double nested loop that uses i for the columns/processes and j for the rows of each column
        if (rank == 0)
        {
            for (int i = 0; i < n_processes; i++)
            {
                for (int j = 0; j < 500; j++)
                {
                    sendColumns[i][j] = values[j][letterColumnMap[argv[i + 3]] - 2];
                }
            }
        }

        // broadcost to processes
        MPI_Bcast(&sendColumns[0][0], 500 * n_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //get operation
        MPIoperation = argv[2];

        if (MPIoperation == "max")
        {
            struct messageData *maxValues = NULL;
            if (rank == 0)
            {
                maxValues = new messageData[n_processes];
            }

            in.value = sendColumns[rank][0];
            for (int j = 1; j < 500; j++)
            {
                if (sendColumns[rank][j] > in.value)
                {
                    in.value = sendColumns[rank][j]; // Store the max value
                    in.index = j; // Store index of max value
                }
            }

            MPI_Gather(&in, 1, MPI_DOUBLE_INT, maxValues, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

            // rank 0 iterates through the results
            if (rank == 0)
            {
                for (int i = 0; i < n_processes; i++)
                {
                    int index = maxValues[i].index;
                    std::cout << MPIoperation << " " << columnHeaders[letterColumnMap[argv[i + 3]]] << " = " << maxValues[i].value << " " << states[index][1] << ", " << states[index][0] << std::endl; // operation name = value state, ST
                }
            }
        }
        else if (MPIoperation == "min")
        {
            struct messageData *minValues = NULL;
            if (rank == 0)
            {
                minValues = (messageData *)malloc(n_processes * sizeof(*minValues));
            }

            in.value = sendColumns[rank][0];
            for (int j = 1; j < 500; j++)
            {
                if (sendColumns[rank][j] < in.value)
                {
                    in.value = sendColumns[rank][j];
                    in.index = j;
                }
            }

            MPI_Gather(&in, 1, MPI_DOUBLE_INT, minValues, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                for (int i = 0; i < n_processes; i++)
                {
                    int index = minValues[i].index;

                    std::cout << MPIoperation << " " << columnHeaders[letterColumnMap[argv[i + 3]]] << " = " << minValues[i].value << " " << states[index][1] << ", " << states[index][0] << std::endl;
                }
            }
        }
        else if (MPIoperation == "avg")
        {
            double avg = 0;
            double *avgA;

            if (rank == 0)
            {
                avgA = new double[n_processes];
            }

            double sum;
            for (int j = 0; j < 500; j++)
            {
                sum += sendColumns[rank][j];
            }
            avg = sum / 500;

            MPI_Gather(&avg, 1, MPI_DOUBLE, avgA, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                for (int i = 0; i < n_processes; i++)
                {
                    std::cout << MPIoperation << " " << columnHeaders[letterColumnMap[argv[i + 3]]] << " = " << avgA[i] << std::endl;
                }
            }
        }
    }

    MPI_Finalize();

    return 0;
}