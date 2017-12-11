#include "graph.h"
#include "tuple.h"
#include <stdio.h>
#include <cli.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdbool.h>


/* Place the code for your Dijkstra implementation in this file */

void
getWeight(component_t c, vertexid_t v1, vertexid_t v2, int* edgeWeight ){
    if (v1 == v2)
        *edgeWeight = 0;


    char s[BUFSIZE];
    ssize_t len, size;
    vertexid_t id1, id2;
    off_t off;
    int readlen, offsetOfInt;
    struct tuple tuple;
    char *buf;
    attribute_t attr;
    bool foundEdge = false;


    /* Open edge file s*/
    memset(s, 0, BUFSIZE);
    sprintf(s, "%s/%d/%d/e", grdbdir, gno, cno);
    c->efd = open(s, O_RDONLY);

    /* Edges */
    size = 0;
    if (c->se != NULL)
        size = schema_size(c->se);

    readlen = (sizeof(vertexid_t) << 1) + size;
    buf = malloc(readlen);

    off = 0;
    while(foundEdge == false){
        lseek(c->efd, off, SEEK_SET);
        len = read(c->efd, buf, readlen);
        if (len <= 0)
            break;

        id1 = *((vertexid_t *) buf);
        id2 = *((vertexid_t *) (buf + sizeof(vertexid_t)));

        if (v1 == id1 && v2 == id2)
            foundEdge = true;

        off += readlen;
    }

    if (foundEdge == false){
        *edgeWeight = -1;
        return;
    }

    memset(&tuple, 0, sizeof(struct tuple));
    tuple.s = c->se;
    tuple.len = size;
    tuple.buf = buf + (sizeof(vertexid_t) << 1);

    for (attr = tuple.s->attrlist; attr != NULL; attr = attr->next) {
        offsetOfInt = tuple_get_offset(&tuple, attr->name);
        if (offsetOfInt >= 0 && attr->bt == INTEGER) {
            *edgeWeight = tuple_get_int(tuple.buf + offsetOfInt);
            break;
        }
        attr = attr->next;
    }

    close(c->efd);
    free(buf);
}


int getNumberOfVertices(){

    char s[BUFSIZE];
    memset(s, 0, BUFSIZE);
    sprintf(s, "%s/%d/%d/v", grdbdir, gno, cno);
    int vfd = open(s, O_RDONLY);

    ssize_t size, len;
    size = sizeof(vertexid_t);

    char* buf = malloc(size);
    off_t off = 0;

    int numberOfVerticies = 0;

    while(true){
        lseek(vfd, off, SEEK_SET);
        len = read(vfd, buf, size);
        if (len <= 0)
            break;

        off += size;
        numberOfVerticies++;
    }

    close(vfd);
    return numberOfVerticies;
}


void getVertices(int numberOfVertices, vertexid_t vertices[numberOfVertices]){

    char s[BUFSIZE];
    memset(s, 0, BUFSIZE);
    sprintf(s, "%s/%d/%d/v", grdbdir, gno, cno);
    int vfd = open(s, O_RDONLY);

    ssize_t size, len;
    size = sizeof(vertexid_t);

    char* buf = malloc(size);
    off_t off = 0;



    for (int i=0; i < numberOfVertices; i++){
        lseek(vfd, off, SEEK_SET);
        len = read(vfd, buf, size);
        if (len <= 0)
            break;

        vertices[i] = *((vertexid_t *) buf);
        off += size;
    }

    close(vfd);
}

int
component_sssp(
        component_t c,
        vertexid_t v1,
        vertexid_t v2,
        int *n,
        int *total_weight,
        vertexid_t **path)
{

    /*
     * Figure out which attribute in the component edges schema you will
     * use for your weight function
     */

    int numberOfVertices = getNumberOfVertices();

    int edgeWeight;

    int D[numberOfVertices];
    int P[numberOfVertices];
    int VminusS[numberOfVertices];
    int indexOfSmallestKnownVertexDistance = -1;
    int indexOfSecondSmallestKnownVertexDistance = -1;
    int indexOfV1 = -1;
    int indexOfV2 = -1;


    vertexid_t V[numberOfVertices];
    getVertices(numberOfVertices, V);

    for (int i = 0; i < numberOfVertices; i++){
        VminusS[i] = 0;
        D[i] = -1;
        if (V[i] == v1){
            indexOfV1 = i;
        }
        if (V[i] == v2){
            indexOfV2 = i;
        }
    }


    D[indexOfV1] = 0;
    int i = indexOfV1;
    while (i != indexOfV2){
        VminusS[i] = -1;
        for (int k = 0; k < numberOfVertices; k++) {
            getWeight(c, V[i], V[k], &edgeWeight);
            if (VminusS[k] != -1 && edgeWeight > -1 && (edgeWeight + D[i] < D[k] || D[k] < 0)) {
                D[k] = edgeWeight + D[i];
                P[k] = i;

                if (indexOfSmallestKnownVertexDistance == -1
                    || (indexOfSmallestKnownVertexDistance > -1 && D[indexOfSmallestKnownVertexDistance] > D[k])) {
                    indexOfSmallestKnownVertexDistance = k;
                } else if (indexOfSecondSmallestKnownVertexDistance == -1
                           || (indexOfSecondSmallestKnownVertexDistance > -1 && D[indexOfSecondSmallestKnownVertexDistance] > D[k])){
                    indexOfSecondSmallestKnownVertexDistance = k;
                }
            }

        }

        if ((indexOfSmallestKnownVertexDistance < indexOfSecondSmallestKnownVertexDistance && indexOfSmallestKnownVertexDistance > -1)
            || indexOfSecondSmallestKnownVertexDistance < 0){
            if (i == indexOfSmallestKnownVertexDistance){
                printf("ERROR: No path found to exist between vertices.\n");
                return -1;
            }

            i = indexOfSmallestKnownVertexDistance;
            indexOfSmallestKnownVertexDistance = -1;
        } else{
            if (i == indexOfSecondSmallestKnownVertexDistance){
                printf("ERROR: No path found to exist between vertices.\n");
                return -1;
            }

            i = indexOfSecondSmallestKnownVertexDistance;
            indexOfSecondSmallestKnownVertexDistance = -1;
        }

    }

    printf("Calculating SSSP\n");

    if (P[indexOfV2] == -1 || D[indexOfV2] == 0){
        printf("ERROR: No path found to exist between vertices.\n");
        return -1;
    }

    *total_weight = D[indexOfV2];
    printf("Total Weight: %d\n", *total_weight);


    int currentVertexIndexOfPath = indexOfV2;
    int numberOfNodes = 1;
    int pathIndexVertices[numberOfVertices];

    for (int i =0; i < numberOfVertices; i++){
        pathIndexVertices[i] = -1;
    }

    P[indexOfV1] = -1;
    pathIndexVertices[0] = indexOfV2;
    while (currentVertexIndexOfPath != indexOfV1){
        currentVertexIndexOfPath = P[currentVertexIndexOfPath];
        pathIndexVertices[numberOfNodes] = currentVertexIndexOfPath;
        numberOfNodes++;
    }

    path[numberOfNodes] = &V[indexOfV1];
    printf("%llu ", *path[numberOfNodes]);
    for (int i = numberOfNodes - 2; i >= 0; i+=-1){
        if (pathIndexVertices[i] == -1) continue;

        path[numberOfNodes] = &V[pathIndexVertices[i]];
        printf("--> %llu ", *path[numberOfNodes]);
    }

    printf("\nNumber of Nodes: %d\n", numberOfNodes);

    *n = numberOfNodes;
    return 0;
}
