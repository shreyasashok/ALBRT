#ifndef NEIGHBOR_RELATIONSHIP_H
#define NEIGHBOR_RELATIONSHIP_H

/**
 * Contains internal information about each octant neighbor.
 * Basically a more minimal version of OctantInternalInfo. 
 * The intended use is to 
*/
struct NeighborInfo {
    int proc;
    int index;
};

enum NeighborType {
    MATCH,
    FINE,
    COARSE,
    BOUNDARY
};

//DEBUG only
static char* neighborTypeStrings[] = {(char*)"MATCH", (char*)"FINE", (char*)"COARSE", (char*)"BOUNDARY"};

/**
 * Describes an octant's relationship to its neighbors. 
 * Uses NeighborInfo to identify neighbor octants.
*/
struct NeighborRelationship {
    NeighborType type;
    union Neighbor {
        NeighborInfo match;
        NeighborInfo fine[4]; //maximum four neighbors for a fine relationship
        struct Coarse {
            NeighborInfo info;
            int index;
        } coarse;
    } neighbor; 
};

#endif