/*
* @file: main.cpp
* @author: Daniel Murga
* @date: 2017.2.20
* @brief: The user interface
*/

#include "Train.h"
#include "barrier.h"
#include <iostream>
#include <fstream>
#include <thread>

int **schedule; // double array that holds the routes for each train at timestep as [Destination][Origin] = 0 or 1
Train *trains; // array that holds train objects created with route, current location, number of stops
std::mutex sMutex; // mutex used to lock schedule and cout
Barrier B; // barrier used to keep threads in sync
char names[26] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'}; // easily name trains by int

/* iterates through the train list to check if all trains have reached destination and returns whether to continue or stop*/
bool end(int totalTrains) {
    for (int i = 0; i < totalTrains; i++)
    {
        if((trains[i].cLocation + 1) != trains[i].nStops) {
            return true;
        }
    }
    return false;
}

/* threads call this function to determine where to send the trains */
void step(int cTrain, int totalTrains, int totalStations) {
    /* make sure all threads start at the same time */
    B.barrier(totalTrains);
    int timestep = 0;
    bool go = true;

    /* iterates until all trains reach destination */
    while(go) { 
        /* check if we're at the end of the route */
        sMutex.lock();
        if((trains[cTrain].cLocation + 1) != trains[cTrain].nStops) {
            /* get the current station and the next station in the route */
            int cStation = trains[cTrain].route[trains[cTrain].cLocation];
            int nStation = trains[cTrain].route[trains[cTrain].cLocation + 1];
            /* check if that track has already been reserved by another train */
            if(schedule[nStation][cStation] == 1 && schedule[cStation][nStation] == 1) {
                std::cout << "At timestep: " << timestep << " train " << names[cTrain] << " must stay at station " << cStation << std::endl;
            }
            /* reserve that track for the train if empty and iterate to next stop */
            else { 
                schedule[nStation][cStation] = 1;
                schedule[cStation][nStation] = 1;
                std::cout << "At timestep: " << timestep << " train " << names[cTrain] << " is going from station " << cStation << " to station " << nStation << std::endl;
                trains[cTrain].cLocation++;
            }
        }
        sMutex.unlock();
        /* make sure all threads have been assigned for that timestep */
        B.barrier(totalTrains);

        /* reset reserved track list */
        if(cTrain == 0) {
            for(int n = 0; n < totalStations; n++) {
                for (int m = 0; m < totalStations; m++)
                {
                    schedule[n][m] = 0;
                }
            }
        }
        timestep++;

        /* check if we should keep going */
        sMutex.lock();
        go = end(totalTrains);
        sMutex.unlock();

        /* make sure the reset did not cause thread 0 to lag behind */
        B.barrier(totalTrains);
    }
}

int main(int argc, char **argv)
{
    int nTrains;
    int nStations;
    int nStops;
    int station;

    /* check for valid arguments */
    if(argc < 2) {
        std::cout << "Pass in a valid argument" << std::endl;
        return(0);
    }

    /* read file and set constants */
    std::ifstream data;
    data.open(argv[1]);
    data >> nTrains;
    data >> nStations;

    /* create arrays to hold threads, trains, and scheduler */
    trains = new Train[nTrains];
    std::thread** threads = new std::thread*[nTrains];
    schedule = new int*[nStations];
    for(int n = 0; n < nStations; n++) {
        schedule[n] = new int[nStations];
        for (int m = 0; m < nStations; m++)
        {
            schedule[n][m] = 0;
        }
    }

    /* fill arrays with data from file */
    std::cout << "Starting simulation..." << std::endl;
    for (int i = 0; i < nTrains; i++)
    {
        /* create route for each train */
        data >> nStops;
        int *route = new int[nStops];
        for (int j = 0; j < nStops; j++)
        {
            data >> station;
            route[j] = station;
        }

        /* create train and place in train array */
        Train T (nStops, route);
        trains[i] = T;

        /* create thread */
        threads[i] = new std::thread(step, i, nTrains, nStations);
    }

    /* wait for all threads to finish and exit program */ 
    for (int l = 0; l < nTrains; l++)
    {
        threads[l]->join();
    }

    return (0);
}