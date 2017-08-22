#ifndef TRAIN_H
#define TRAIN_H

class Train
{
  public:
    int cLocation;
    int nStops;
    int *route;

    Train(){}
    Train(int stops, int *rt) {
        cLocation = 0;
        nStops = stops;
        route = rt;
    }
};

#endif