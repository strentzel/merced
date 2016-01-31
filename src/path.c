#include <assert.h>
#include <stdlib.h>

#include "dscan.h"
#include "dscanimp.h"

static vector2_t compute_normal( vector2_t s, vector2_t t )
{
    vector2_t v = minus2( t, s );
    vector2_t r;
    r.x = -v.y;
    r.y =  v.x;
    return normalize2( r );
}

static vector2_t get_init_point( convexpolygon_t* poly, projection_t* proj ) {
    vector2_t work = scale2( poly->mVecs[proj->edge_start],proj->scale_start );
    work = add2( poly->mPoints[proj->edge_start], work );
    return work;
}

static vector2_t get_final_point( convexpolygon_t* poly, projection_t* proj )
{
    vector2_t work = scale2( poly->mVecs[proj->edge_end], proj->scale_end );
    work = add2( poly->mPoints[proj->edge_end], work);
    return work;
}

static projection_t project( convexpolygon_t* poly, 
                             vector2_t surface_start, 
                             vector2_t surface_end,
                             vector2_t ray,
                             lightpath_type_t type )
{
    if( LIGHTPATH_EXTERNAL == type )
    {
        return convexpolygon_project( poly, surface_start, surface_end, ray );
    }
    else
    {
        assert( LIGHTPATH_INTERNAL == type );
        return convexpolygon_project( poly, surface_end, surface_start, 
                                      scale2( ray, -1.0F ) );
    }
}

static vector2_t* get_points( convexpolygon_t* poly, projection_t* proj, 
                              lightpath_type_t type, int* cnt )
{
    int n =  (proj->edge_end-proj->edge_start+poly->mPointsCnt) %
        poly->mPointsCnt + 2;
    vector2_t* retval = malloc( sizeof(vector2_t) * n );
    int i;
    int p;

    if( proj->edge_start == -1 )
    {
        *cnt = 0;
        return retval;
    }

    retval[0] = get_init_point( poly, proj );
    retval[n-1] = get_final_point( poly, proj );
    p = (proj->edge_start + 1) % poly->mPointsCnt;
    for( i = 0; i < n-2; ++i, p = (p+1)%poly->mPointsCnt )
    {
        retval[i+1] = poly->mPoints[p];
    }
    
    *cnt = n;
    return retval;
}

static vector2_t reverse_project( vector2_t point, vector2_t start,
                                  vector2_t end, vector2_t ray )
{
    vector2_t rel = minus2( point, start );
    vector2_t vec = normalize2( minus2( end, start ) );
    return add2( start, scale2( vec, dot2( vec, rel ) ) );
}

static lightpath_node_t* make_surface_node( vector2_t start, vector2_t end,
                                            vector2_t ray, float weight,
                                            lightpath_node_t* next )
{
    lightpath_node_t* retval =  malloc( sizeof(lightpath_node_t) );
    retval->next = next;
    retval->data.surface = malloc( sizeof(lightpath_surface_t) );
    retval->data.surface->surface_start = start;
    retval->data.surface->surface_end = end;
    retval->data.surface->ray = ray;
    retval->weight = weight;
    retval->type = LIGHTPATH_INTERNAL;
    return retval;
}

static lightpath_node_t* make_leaf_node( vector2_t start, vector2_t end,
                                         vector2_t ray, float weight,
                                            lightpath_node_t* next )
{
    lightpath_node_t* retval =  malloc( sizeof(lightpath_node_t) );
    retval->next = next;
    retval->data.surface = malloc( sizeof(lightpath_surface_t) );
    retval->data.surface->surface_start = start;
    retval->data.surface->surface_end = end;
    retval->data.surface->ray = ray;
    retval->weight = weight;
    retval->type = LIGHTPATH_LEAF;
    return retval;
}

static lightpath_polygon_t* make_shape( vector2_t p0, vector2_t p1, 
                                        vector2_t* p, int pn )
{
    lightpath_polygon_t* retval;
    int i;
    retval = malloc( sizeof(lightpath_polygon_t) );
    retval->points = malloc( sizeof(vector2_t) * (pn+2) );
    retval->points[0] = p0;
    retval->points[1] = p1;
    for( i = 0; i < pn; ++i )
    {
        retval->points[i+2] = p[i];
    }
    retval->pointsCnt = pn+2;
    return retval;
}
		
static void estimate_refraction( material_t* outside, material_t* inside,
                                 vector2_t start, vector2_t end,
                                 vector2_t incoming, 
                                 vector2_t* refr, float* refr_weight,
                                 float* density_weight )
{
    vector2_t surfaceNormal = compute_normal( start, end );
    float nOut = (outside->n_max + outside->n_min) / 2.0F;
    float sinThetaOut = min( 1.0F, max( -1.0F, 
                                        -cross2( surfaceNormal, incoming ) / 
                                        (norm2( surfaceNormal ) * 
                                         norm2( incoming )) ) );
    float nIn = (inside->n_max + inside->n_min) / 2.0F;
    float sinThetaIn = min( 1.0F, max( -1.0F, sinThetaOut * nOut / nIn ) );
    float cosThetaIn = sqrt(1 - sinThetaIn*sinThetaIn);
    float cosThetaOut = sqrt( 1 - sinThetaOut*sinThetaOut );
    refr->x = cosThetaIn*(-surfaceNormal.x) + -sinThetaIn*(-surfaceNormal.y);
    refr->y = sinThetaIn*(-surfaceNormal.x) + cosThetaIn*(-surfaceNormal.y); 
		
    if( cosThetaOut < outside->cos )
    {
        *refr_weight = 0.0F;
        *density_weight = 0.0F;
    }
    else
    {
        float nIn = (inside->n_max + inside->n_min) / 2.0F;
        float s = (nOut*cosThetaOut - nIn*cosThetaIn)/
            (nOut*cosThetaOut + nIn*cosThetaIn);
        float Rs = s * s;
        float p = (nOut*cosThetaIn - nIn*cosThetaOut)/
            (nOut*cosThetaIn + nIn*cosThetaOut);
        float Rp = p * p;
        float T = 1.0F - (Rs + Rp) / 2.0F;
        float maxTd = 2 * nOut + 2 / nOut;
        *refr_weight = T;
        if( T * fabs(cosThetaOut) < maxTd * fabs(cosThetaIn) )
        {
            *density_weight = T * cosThetaOut / cosThetaIn;
        }
        else
        {
            *density_weight = maxTd;
        }
    }
}

static vector2_t compute_reflection( vector2_t start, vector2_t end,
                                     vector2_t ray )
{
    vector2_t surface_normal = compute_normal( start, end );
    return add2( scale2( surface_normal, -2.0F * dot2( surface_normal, ray ) ),
                 ray );
}

lightpath_t* dscan_lightpath_create
(
 vector2_t source_start,
 vector2_t source_end
)
{
    lightpath_t* retval = malloc( sizeof(lightpath_t) );
    vector2_t dir = minus2( source_end, source_start );
    retval->head = retval->tail = malloc( sizeof(lightpath_row_t) );
    retval->head->head = malloc( sizeof(lightpath_node_t) );
    retval->head->head->data.surface = malloc( sizeof(lightpath_surface_t) );
    retval->head->head->data.surface->surface_start = source_start;
    retval->head->head->data.surface->surface_end = source_end;
    retval->head->head->data.surface->ray.x = -dir.y;
    retval->head->head->data.surface->ray.y =  dir.x;
    retval->head->head->weight = 1.0F;
    retval->head->head->type = LIGHTPATH_EXTERNAL;
    retval->head->head->next = NULL;
    retval->head->next = NULL;
    retval->point_cnt = 0;
    retval->element_cnt = 0;
    return retval;
} 

int dscan_convexpolygon_trace
(
 lightpath_t* path,
 material_t* external,
 material_t* internal,
 convexpolygon_t* cross_section,
 float minimum_weight,
 float minimum_surface
)
{
    lightpath_row_t* row = path->tail;
    lightpath_row_t* next;
    lightpath_node_t* node;
    lightpath_node_t* leafs = NULL;
    bool cont = 0;

    if( NULL == row ) return 0;

    next = malloc( sizeof(lightpath_row_t) );
    next->head = NULL;

    node = row->head;
    assert( node != NULL );
    while(1)
    {
        projection_t segment;
        vector2_t start = node->data.surface->surface_start;
        vector2_t end = node->data.surface->surface_end;
        vector2_t* pp;
        int ppCnt;
    
        assert( LIGHTPATH_EXTERNAL == node->type ||
                LIGHTPATH_INTERNAL == node->type );

        segment = project( cross_section, start, end, 
                           node->data.surface->ray, node->type );
        pp = get_points( cross_section, &segment, node->type, &ppCnt );

        if( LIGHTPATH_EXTERNAL == node->type && ppCnt > 1 )
        {
            vector2_t startP;
            int i;

            startP = reverse_project( pp[0], start, end, 
                                      node->data.surface->ray );
            end = reverse_project( pp[ppCnt-1], start, end,
                                   node->data.surface->ray);
            start = startP;

            for( i = 0; i < ppCnt-1; ++i )
            {
                vector2_t refr;
                float refr_weight;
                float unused;
                
                if( norm2( minus2( pp[i+1], pp[i] ) ) < minimum_surface )
                    continue;

                estimate_refraction( external, internal,
                                     pp[i+1], pp[i], node->data.surface->ray,
                                     &refr, &refr_weight, &unused );
                if( node->weight * refr_weight >= minimum_weight )
                {
                    next->head = make_surface_node( pp[i], pp[i+1], refr,
                                                    node->weight * refr_weight,
                                                    next->head );
                    cont = 1;
                }
            }
        }
        else
        {
            int i;

            for( i = 0; i < ppCnt-1; ++i )
            {
                vector2_t refr, refl;
                float refr_weight, refl_weight, density_weight;

                if( norm2( minus2( pp[i+1], pp[i] ) ) < minimum_surface )
                    continue;

                estimate_refraction( internal, external,
                                     pp[i], pp[i+1], node->data.surface->ray,
                                     &refr, &refr_weight, &density_weight );
                refl = compute_reflection( pp[i], pp[i+1], 
                                           node->data.surface->ray);
                refl_weight = 1 - refr_weight;

                if( node->weight * refl_weight >= minimum_weight )
                {
                    next->head = make_surface_node( pp[i], pp[i+1], refl,
                                                    node->weight*refl_weight,
                                                    next->head );
                    cont = 1;
                }
                if( node->weight * density_weight >= minimum_weight )
                {
                    leafs = make_leaf_node( pp[i+1], pp[i], 
                                            node->data.surface->ray,
                                            node->weight*refr_weight, leafs );
                    path->element_cnt++;
                    path->point_cnt += 3;
                }
            }
        }

        free( node->data.surface );
        node->data.polygon = make_shape( start, end, pp, ppCnt );
        node->type = LIGHTPATH_POLYGON;
        path->element_cnt++;
        path->point_cnt += node->data.polygon->pointsCnt;
        free( pp );

        if( node->next != NULL )
            node = node->next;
        else
            break;
    }
    
    node->next = leafs;

    if( cont )
    {
        row->next = next;
        path->tail = next;
        next->next = NULL;
        return 1;
    }
    else
    {
        free( next );
        return 0;
    }
}
    
int dscan_lightpath_points(lightpath_t* lp)
{
    return lp->point_cnt;
}

int dscan_lightpath_elements(lightpath_t* lp)
{
    return lp->element_cnt;
}

void dscan_lightpath_marshall(lightpath_t* lp, 
                              int element_first, int element_last,
                              dscan_vector2_t* points, int* element_pointCnts, 
                              int* element_types, float* element_weights)
{
    int cnt = 0;
    lightpath_row_t* row = lp->head;
    while( row != NULL )
    {
        lightpath_node_t* node = row->head;
        while( node != NULL )
        {
            if( LIGHTPATH_POLYGON == node->type )
            {
                if( cnt >= element_first && cnt < element_last )
                {
                    int j;
                    for( j = 0; j < node->data.polygon->pointsCnt; ++j )
                    {
                        *(points++) = node->data.polygon->points[j];
                    }
                    *(element_pointCnts++) = node->data.polygon->pointsCnt;
                    *(element_types++) = 0;
                    *(element_weights++) = node->weight;
                }
                ++cnt;
            }
            else
            {
                if( cnt >= element_first && cnt < element_last )
                {
                    assert( LIGHTPATH_LEAF == node->type );
                    *(points++) = node->data.surface->surface_start;
                    *(points++) = node->data.surface->surface_end;
                    *(points++) = node->data.surface->ray;
                    *(element_pointCnts++) = 3;
                    *(element_types++) = 1;
                    *(element_weights++) = node->weight;
                }
                ++cnt;
            }
            node = node->next;
        }
        row = row->next;
    }
}

void dscan_lightpath_destroy
(
 lightpath_t* path
)
{
    lightpath_row_t* row = path->head;
    while( row != NULL )
    {
        lightpath_row_t* next = row->next;
        lightpath_node_t* node = row->head;
        while( node != NULL )
        {
            lightpath_node_t* next = node->next;
            if( LIGHTPATH_POLYGON == node->type )
            {
                free( node->data.polygon->points );
                free( node->data.polygon );
            }
            else
            {
                free( node->data.surface );
            }
            free( node );
            node = next;
        }
        free( row );
        row = next;
    }
    free( path );
}
