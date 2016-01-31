#ifndef DSCANIMP_H__
#define DSCANIMP_H__

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265F
#endif

#ifndef bool
#define bool char
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef min
#define min(a,b) (a < b ? a : b)
#endif

#ifndef max
#define max(a,b) (a < b ? b : a)
#endif

typedef struct dscan_vector2 vector2_t;
typedef struct dscan_matrix3 matrix3_t;
typedef struct dscan_roundcut_measurements roundcut_measurements_t;
typedef struct dscan_roundcut_errors roundcut_errors_t;
typedef struct dscan_material material_t;
typedef struct dscan_observer observer_t;
typedef struct dscan_lightpath lightpath_t;

struct dscan_vector3
{
    float x;
    float y;
    float z;
};
typedef struct dscan_vector3 vector3_t;

struct dscan_convexpolygon
{
    vector2_t* mPoints;
    vector2_t* mVecs;
    int mPointsCnt;
    vector2_t center;
    float radius2;
};
typedef struct dscan_convexpolygon convexpolygon_t;

struct dscan_edge
{
    int a;
    int b;
};
typedef struct dscan_edge edge_t;

struct dscan_face
{
    int* edges;
    int edgesCnt;
    char reverse0;
    char reverse1;
};
typedef struct dscan_face face_t;

struct dscan_convexpolyhedron
{
    vector3_t* mOrigPoints;
    vector3_t* mPoints;
    int mPointsCnt;
    edge_t* mEdges;
    int mEdgesCnt;
    face_t* mFaces;
    int mFacesCnt;
};
typedef struct dscan_convexpolyhedron convexpolyhedron_t;

struct dscan_projection
{
    int edge_start;
    float scale_start;
    int edge_end;
    float scale_end;
    vector2_t source_start;
    vector2_t source_vec;
};
typedef struct dscan_projection projection_t;

struct dscan_lightpath_surface
{
    vector2_t surface_start;
    vector2_t surface_end;
    vector2_t ray;    
};
typedef struct dscan_lightpath_surface lightpath_surface_t;

struct dscan_lightpath_polygon
{
    vector2_t* points;
    int pointsCnt;
};
typedef struct dscan_lightpath_polygon lightpath_polygon_t;

enum dscan_lightpath_type
{
    LIGHTPATH_EXTERNAL,
    LIGHTPATH_INTERNAL,
    LIGHTPATH_POLYGON,
    LIGHTPATH_LEAF
};
typedef enum dscan_lightpath_type lightpath_type_t;

struct dscan_lightpath_node
{
    union
    {
        lightpath_polygon_t* polygon;
        lightpath_surface_t* surface;
    } data;
    float weight;
    lightpath_type_t type;
    struct dscan_lightpath_node* next;
};
typedef struct dscan_lightpath_node lightpath_node_t;

struct dscan_lightpath_row
{
    lightpath_node_t* head;
    struct dscan_lightpath_row* next;
};
typedef struct dscan_lightpath_row lightpath_row_t;

struct dscan_lightpath
{
    lightpath_row_t* head;
    lightpath_row_t* tail;
    int point_cnt;
    int element_cnt;
};

struct dscan_material
{
    float n_min;
    float n_max;
    float theta;
    float cos;

    float c0;
    float c1;
    float c2;
};

struct dscan_observer
{
    int type;
};

struct dscan_prism
{
    material_t* material_incident;
    material_t* material_transmittent;
    observer_t* observer;
    float light_weight;

    vector2_t mSurfacePoint;
    vector2_t mSurfaceVec;
    vector2_t mMinDir;
    vector2_t mMaxDir;

    float mThetaD, mSinThetaD, mCosThetaD;
    float mSinMinTheta, mSinMaxTheta;

    vector2_t mNegB;
    vector2_t mNegC;

    vector2_t mNormA;
    vector2_t mNormB;
    vector2_t mNormC;
    vector2_t mNormD;
    vector2_t mNormE;

    bool mHasCoreSample;
    int mCoreRgb;
};
typedef struct dscan_prism prism_t;

vector3_t normalize3(vector3_t v);
vector3_t add3(vector3_t a, vector3_t b);
vector3_t minus3(vector3_t a, vector3_t b);
vector3_t scale3(vector3_t v, float s);
vector3_t cross3(vector3_t a, vector3_t b);
vector3_t intersect_plane(vector3_t norm, vector3_t rpoint, vector3_t rvec);
float dot3(vector3_t s, vector3_t t);
float norm3(vector3_t a);
vector2_t multiply3(matrix3_t m, vector2_t v, bool z);
vector2_t normalize2(vector2_t v);
vector2_t add2(vector2_t a, vector2_t b);
vector2_t minus2(vector2_t a, vector2_t b);
vector2_t scale2(vector2_t v, float s);
float cross2(vector2_t a, vector2_t b);
float dot2(vector2_t a, vector2_t b);
float norm2(vector2_t a);
float compute_intersect_scale( vector2_t segment, vector2_t rel_point, vector2_t rel_vec, bool* err );
projection_t convexpolygon_project( convexpolygon_t* poly,
                                    vector2_t surface_start,
                                    vector2_t surface_end,
                                    vector2_t ray );
convexpolygon_t convexpolygon_intersect( convexpolygon_t* s, convexpolygon_t* t );
convexpolygon_t convexpolygon_init( vector2_t* points, int cnt );
void convexpolygon_destroy( convexpolygon_t poly );
void convexpolygon_check( convexpolygon_t* poly );
convexpolyhedron_t convexpolyhedron_init(vector3_t* points, int pointsCnt,
                                         edge_t* edges, int edgesCnt,
                                         face_t* faces, int facesCnt);
void convexpolyhedron_check(convexpolyhedron_t* poly, bool verbose);


#endif /*DSCANIMP_H__*/
