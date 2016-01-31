#ifndef DSCAN_H__
#define DSCAN_H__

struct dscan_vector2
{
    float x;
    float y;
};
typedef struct dscan_vector2 dscan_vector2_t;

struct dscan_matrix3
{
    float d[9];
};
typedef struct dscan_matrix3 dscan_matrix3_t;

struct dscan_roundcut_input
{
    float* width_min;
    float* width_max;
    float* height;
    float* table_pct;
    float* crownheight_pct;
    float* pavilionheight_pct;
    float* pavilion_angle;
    float* width_avg;
    float* height_pct;
    float* crown_angle;
    float* girdle_pct;
    float* culet_pct;
    float* star_pct;
    float* pavilion_pct;
    int* girdleedges;
};
typedef struct dscan_roundcut_input dscan_roundcut_input_t;

struct dscan_roundcut_measurements
{
    float table_pct;
    float crownheight_pct;
    float pavilionheight_pct;
    float girdle_pct;
    float star_pct;
    float pavilion_pct;
    float culet_pct;
    int girdleedges;
};
typedef struct dscan_roundcut_measurements dscan_roundcut_measurements_t;

struct dscan_roundcut_errors
{
    float min_girdle_pct;
};
typedef struct dscan_roundcut_errors dscan_roundcut_errors_t;

typedef struct dscan_convexpolygon* dscan_convexpolygon_t;
typedef struct dscan_convexpolyhedron* dscan_convexpolyhedron_t;
typedef struct dscan_material* dscan_material_t;
typedef struct dscan_observer* dscan_observer_t;
typedef struct dscan_lightpath* dscan_lightpath_t;
typedef struct dscan_prism* dscan_prism_t;

int dscan_roundcut_measure(dscan_roundcut_input_t in, dscan_roundcut_input_t* out);
dscan_roundcut_measurements_t dscan_roundcut_complete(dscan_roundcut_input_t given, dscan_roundcut_input_t calculated);
dscan_convexpolyhedron_t dscan_roundcut_create(dscan_roundcut_measurements_t measurements, dscan_roundcut_errors_t* err);
int dscan_roundcut_rotations();
dscan_matrix3_t dscan_roundcut_rotation(int index);
dscan_convexpolygon_t dscan_roundcut_crosssection(dscan_convexpolyhedron_t poly, int index);

dscan_material_t dscan_material_create(const char* type, int len);
void dscan_material_destroy(dscan_material_t material);

dscan_observer_t dscan_observer_create(const char* type, int len);
void dscan_observer_destroy(dscan_observer_t observer);

void dscan_convexpolyhedron_destroy(dscan_convexpolyhedron_t poly);
void dscan_convexpolyhedron_transform(dscan_convexpolyhedron_t poly, dscan_matrix3_t transform);
int dscan_convexpolyhedron_edges(dscan_convexpolyhedron_t poly);
float dscan_convexpolyhedron_volume(dscan_convexpolyhedron_t poly);
int dscan_convexpolyhedron_visible(dscan_convexpolyhedron_t poly, float* edge_points);

dscan_convexpolygon_t dscan_convexpolygon_create(dscan_vector2_t* points, int pointsCnt);
void dscan_convexpolygon_destroy(dscan_convexpolygon_t poly);
int dscan_convexpolygon_points(dscan_convexpolygon_t poly);
void dscan_convexpolygon_marshall(dscan_convexpolygon_t poly, float* points);
int dscan_convexpolygon_trace(dscan_lightpath_t path, dscan_material_t external, dscan_material_t internal, dscan_convexpolygon_t cross_section, float minimum_weight, float minimum_surface );

dscan_lightpath_t dscan_lightpath_create(dscan_vector2_t source_start, dscan_vector2_t source_end);
void dscan_lightpath_destroy(dscan_lightpath_t lp);
int dscan_lightpath_points(dscan_lightpath_t lp);
int dscan_lightpath_elements(dscan_lightpath_t lp);
void dscan_lightpath_marshall(dscan_lightpath_t lp, int element_first, int element_last, dscan_vector2_t* points, int* element_pointCnts, int* element_types, float* element_weights );

dscan_prism_t dscan_prism_create( dscan_vector2_t source_start, dscan_vector2_t source_end, dscan_vector2_t incidence, dscan_material_t material_incident, dscan_material_t material_transmitent, dscan_observer_t observer, dscan_matrix3_t display_transform, float weight );
void dscan_prism_destroy( dscan_prism_t prism );
void dscan_prism_rasterize( dscan_prism_t prism, int* bitmap, int offset, int stride, int w, int h);

#endif /*DSCAN_H__*/
