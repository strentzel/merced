#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dscan.h"
#include "dscanimp.h"

#ifdef CATCHEXIT
void catch_exit(int);
#define EXIT catch_exit
#else
#define EXIT exit
#endif

static float ComputeIntersectScale( vector2_t s, vector2_t rp, vector2_t v )
{
    bool err;
    return compute_intersect_scale( s, rp, v, &err );
}

static convexpolygon_t Empty()
{
    convexpolygon_t retval;

    retval.mPoints = 0;
    retval.mVecs = 0;
    retval.mPointsCnt = 0;
    retval.center.x = 0;
    retval.center.y = 0;
    retval.radius2 = 0;

    return retval;
}

static void GetReverse(convexpolyhedron_t* cp, bool* reverse, int cnt, int face) {
    int i;

    assert(cnt >= cp->mFaces[face].edgesCnt);
    reverse[0] = cp->mFaces[face].reverse0;
    reverse[1] = cp->mFaces[face].reverse1;
    for( i = 2; i < cp->mFaces[face].edgesCnt; ++i) {
        int start = (reverse[i - 1] 
                     ? cp->mEdges[cp->mFaces[face].edges[i - 1]].a 
                     : cp->mEdges[cp->mFaces[face].edges[i - 1]].b);
        reverse[i] = (start != cp->mEdges[cp->mFaces[face].edges[i]].a);
    }
}

static bool* MakeReverse(convexpolyhedron_t* cp, int face) 
{
    bool* retval = malloc(sizeof(bool)*cp->mFaces[face].edgesCnt);
    GetReverse(cp, retval, cp->mFaces[face].edgesCnt, face);
    return retval;
}

static int GetVertexIndexR(convexpolyhedron_t* cp, int face, int point, bool reverse) 
{
    int idx = reverse 
        ? cp->mEdges[cp->mFaces[face].edges[point]].b 
        : cp->mEdges[cp->mFaces[face].edges[point]].a;
    return idx;
}

static vector3_t GetVertexR(convexpolyhedron_t* cp, int face, int point, bool reverse) 
{
    return cp->mPoints[GetVertexIndexR( cp, face, point, reverse )];
}

static int GetVertexIndex(convexpolyhedron_t* cp, int face, int point)
{
    bool* reverse = MakeReverse(cp, face);
    int retval = GetVertexIndexR(cp, face, point, reverse[point]);
    free(reverse);
    return retval;
}

static vector3_t GetVertex(convexpolyhedron_t* cp, int face, int point)
{
    return cp->mPoints[GetVertexIndex(cp, face, point)];
}

static vector3_t makeNormal(vector3_t a, vector3_t b, vector3_t c)
{
    float x0 = b.x - a.x;
    float y0 = b.y - a.y;
    float z0 = b.z - a.z;
    float x1 = c.x - b.x;
    float y1 = c.y - b.y;
    float z1 = c.z - b.z;
    vector3_t result;
    result.x = y0 * z1 - y1 * z0;
    result.y = z0 * x1 - z1 * x0;
    result.z = x0 * y1 - x1 * y0;
    return normalize3( result );
}

static vector3_t GetNorm(convexpolyhedron_t* cp, int face)
{
    vector3_t a = GetVertex(cp, face, 0);
    vector3_t b = GetVertex(cp, face, 1);
    vector3_t c = GetVertex(cp, face, 2);
    vector3_t normal = makeNormal(a, b, c);
    return normal;
}

static bool IsVisible(convexpolyhedron_t* poly, int f)
{
    vector3_t p0 = poly->mFaces[f].reverse0 
        ? poly->mPoints[poly->mEdges[poly->mFaces[f].edges[0]].b] 
        : poly->mPoints[poly->mEdges[poly->mFaces[f].edges[0]].a];
    vector3_t p1 = poly->mFaces[f].reverse0 
        ? poly->mPoints[poly->mEdges[poly->mFaces[f].edges[0]].a] 
        : poly->mPoints[poly->mEdges[poly->mFaces[f].edges[0]].b];
    vector3_t p2 = poly->mFaces[f].reverse1 
        ? poly->mPoints[poly->mEdges[poly->mFaces[f].edges[1]].a] 
        : poly->mPoints[poly->mEdges[poly->mFaces[f].edges[1]].b];
    float px = p1.x - p0.x;
    float py = p1.y - p0.y;
    float qx = p2.x - p1.x;
    float qy = p2.y - p1.y;
    return px * qy - py * qx > 0;
}

convexpolygon_t convexpolygon_init
(
 vector2_t* points,
 int pointsCnt
)
{
    convexpolygon_t retval;
    int i;
    float xsum, ysum;
    float r2;

    if( 0 == pointsCnt )
    {
        assert( 0 == points );
        return Empty();
    }

    retval.mPoints = points;
    retval.mPointsCnt = pointsCnt;
    retval.mVecs = malloc( sizeof(vector2_t) * pointsCnt );
    for( i = 0; i < pointsCnt; ++i )
    {
        int next = (i + 1) % pointsCnt;
        retval.mVecs[i].x = points[next].x - points[i].x;
        retval.mVecs[i].y = points[next].y - points[i].y;
    }

    xsum = 0.0F;
    ysum = 0.0F;
    for( i = 0; i < pointsCnt; ++i ) {
        xsum += points[i].x;
        ysum += points[i].y;
    }
    retval.center.x = xsum / pointsCnt;
    retval.center.y = ysum / pointsCnt;

    r2 = 0.0F;
    for( i = 0; i < pointsCnt; ++i )
    {
        float x = retval.mPoints[i].x - retval.center.x;
        float y = retval.mPoints[i].y - retval.center.y;
        float c2 = x*x + y*y;
        if( c2 > r2 )
        {
            r2 = c2;
        }
    }
    retval.radius2 = r2;

    return retval;
}

convexpolygon_t* dscan_convexpolygon_create( vector2_t* points, int cnt )
{
    convexpolygon_t* retval = malloc( sizeof(convexpolygon_t) );
    *retval = convexpolygon_init( points, cnt );
    return retval;
}

void dscan_convexpolygon_destroy
(
 convexpolygon_t* poly
)
{
    free( poly->mPoints );
    free( poly->mVecs );
    free( poly );
}

void convexpolygon_destroy
(
 convexpolygon_t poly
)
{
    free( poly.mPoints );
    free( poly.mVecs );
}

int dscan_convexpolygon_points
(
 dscan_convexpolygon_t poly
)
{
    return poly->mPointsCnt;
}

void dscan_convexpolygon_marshall
(
 dscan_convexpolygon_t poly,
 float* points
)
{
    int i;
    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        points[2*i] = poly->mPoints[i].x;
        points[2*i+1] = poly->mPoints[i].y;
    }
}

static void fail2()
{
    fprintf( stderr, "Hard error in convexpolygon_check.\n" );
    EXIT(-1);
}

void convexpolygon_check
(
 convexpolygon_t* poly
)
{
    int i, j;

    if( 0 == poly->mPointsCnt )
    {
        return;
    }

    if( poly->mPointsCnt <= 2 )
    {
        fprintf( stderr, "Too few points.\n" );
        fail2();
    }

    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        int next = (i + 1) % poly->mPointsCnt;
        if( cross2( poly->mVecs[i], poly->mVecs[next] ) <= 0 )
        {
            fprintf( stderr, "Non-convex.\n" );
            fail2();
        }

        for( j = 0; j < poly->mPointsCnt; ++j )
        {
            if( i != j )
            {
                vector2_t r;
                r.x = poly->mPoints[j].x - poly->mPoints[i].x;
                r.y = poly->mPoints[j].y - poly->mPoints[i].y;
                if( cross2( poly->mVecs[i], r ) < 0 )
                {
                    fprintf( stderr, "Not simple.\n" );
                    fail2();
                }
            }
        }
    }
}

static projection_t projection_empty()
{
    projection_t retval;
    retval.edge_start = -1;
    retval.scale_start = 0;
    retval.edge_end = -1;
    retval.scale_end = 0;
    retval.source_start.x = 0;
    retval.source_start.y = 0;
    retval.source_vec.x = 0;
    retval.source_vec.y = 0;
    return retval;
}
	
static projection_t projection_init
(
 int es,
 float ss, 
 int ee,
 float se,
 vector2_t ps,
 vector2_t pv
)
{
    projection_t retval;
    retval.edge_start = es;
    retval.scale_start = ss;
    retval.edge_end = ee;
    retval.scale_end = se;
    retval.source_start = ps;
    retval.source_vec = normalize2( pv );
    return retval;
}

bool is_outside
(
 vector2_t point,
 vector2_t point_on_line,
 vector2_t ray
)
{
    vector2_t rel = minus2( point, point_on_line );
    return cross2( rel, ray ) > 0;
}

projection_t convexpolygon_project
(
 convexpolygon_t* poly,
 vector2_t surface_start,
 vector2_t surface_end,
 vector2_t ray
)
{
    int rear_point;
    float rear_dist;
    int edge0, edge1;
    float scale0, scale1;
    bool found;
    int i;

    if( 0 == poly->mPointsCnt ) return projection_empty();

    rear_point = 0;
    rear_dist = dot2( poly->mPoints[0], ray );
    for( i = 1; i < poly->mPointsCnt; ++i )
    {
        float dist = dot2( poly->mPoints[i], ray );
        if( dist > rear_dist )
        {
            rear_point = i;
            rear_dist = dist;
        }
    }

    /* Find Start */
    edge0 = rear_point;
    found = 0;
    do
    {
        int next = (edge0 + 1) % poly->mPointsCnt;
        if( cross2( poly->mVecs[edge0], ray ) > 0 &&
            is_outside( poly->mPoints[next], surface_start, ray ) )
        {
            found = 1;
            break;
        }
        
        edge0 = next;
    } while( edge0 != rear_point );
    
    if( !found ||
        is_outside( poly->mPoints[edge0], surface_end, ray ) )
    {
        return projection_empty();
    }

    /* Find End */
    edge1 = edge0;
    do
    {
        if( cross2( poly->mVecs[edge1], ray ) <= 0 ||
            is_outside( poly->mPoints[edge1], surface_end, ray ) )
        {
            break;
        }
        
        edge1 = (edge1 + 1) % poly->mPointsCnt;
    } while( edge1 != rear_point );

    edge1 = (edge1 + poly->mPointsCnt - 1) % poly->mPointsCnt;
    
    scale0 = ComputeIntersectScale( poly->mVecs[edge0],
              minus2( surface_start, poly->mPoints[edge0] ), ray );
    scale1 = ComputeIntersectScale( poly->mVecs[edge1],
              minus2( surface_end, poly->mPoints[edge1] ), ray );
    assert( scale0 >= 0 && scale0 <= 1 );
    assert( scale1 >= 0 && scale1 <= 1 );
    if( edge0 == edge1 && scale1 <= scale0 ) 
        return projection_empty();
    else 
        return projection_init( edge0, scale0,
                                edge1, scale1,
                                surface_start,
                                minus2( surface_end, surface_start ) );                                
}

float DistanceFromLine( vector2_t p, vector2_t orig, vector2_t dir )
{
    vector2_t b = minus2( p, orig );
    vector2_t a = normalize2( dir );
    vector2_t x = minus2( b, scale2( a, dot2( a, b ) ) );
    return norm2( x );
}

float ComputeSide( vector2_t p, vector2_t o, vector2_t v )
{
    vector2_t r = minus2( p, o );
    return cross2( v, r );    
}

vector2_t ComputeHalfPlaneIntersect( vector2_t p, vector2_t v,
                                     vector2_t orig, vector2_t dir,
                                     bool* err )
{
    vector2_t rp, retval;
    float scale;
    
    rp = minus2( orig, p );
    scale = compute_intersect_scale( v, rp, dir, err );
    retval = add2( p, scale2( v, scale ) );
    return retval;
}

void Rotate( vector2_t* p, int cnt, int f )
{
    int first, middle, next;
    vector2_t tmp;

    first = 0;
    middle = f;
    next = f;
    while( first != next )
    {
        tmp = p[first];
        p[first++] = p[next];
        p[next++] = tmp;
        if( next == cnt ) next = middle;
        if( first == middle ) middle = next;
    }
}

void IntersectHalfPlane
(
 vector2_t* points,
 int* cnt,
 vector2_t orig,
 vector2_t dir
)
{
    int i, f, l;
    vector2_t fv, lv;
    vector2_t fp, lp;
    bool ferr = 0, lerr = 0;
    float mx = 0, mn = 0;
    int mxi = -1, mni = -1;
    
    for( i = 0; i < *cnt; ++i )
    {
        float s = ComputeSide( points[i], orig, dir );
        float d = DistanceFromLine( points[i], orig, dir ) *
            (s > 0 ? 1 : (s < 0 ? -1 : 0));
        if( d > mx )
        {
            mx = d;
            mxi = i;
        }
        if( d < mn )
        {
            mn = d;
            mni = i;
        }
    }

    if( mx <= 0 )
    {
        *cnt = 0;
        return;
    }
    assert( !(mn < 0) || ComputeSide( points[mni], orig, dir ) < 0 );
    if( mn >= 0 )
    {
        return;
    }

    for( l = mxi; ComputeSide(points[l],orig,dir)>=0; l = (l+1)%*cnt );
    for( f = mxi; ComputeSide(points[f],orig,dir)>=0; f = (f+*cnt-1)%*cnt );
    f = (f+1)%*cnt;

    fv = minus2( points[(f+*cnt-1)%*cnt], points[f] );
    fp = ComputeHalfPlaneIntersect( points[f], fv, orig, dir, &ferr );
    lv = minus2( points[l], points[(l+*cnt-1)%*cnt] );
    lp = ComputeHalfPlaneIntersect( points[(l+*cnt-1)%*cnt], lv, orig, dir, &lerr );
    
    Rotate( points, *cnt, f );
    *cnt = (f<l ? l-f : l+*cnt - f);
    if( !lerr )
    {
        points[(*cnt)++] = lp;
    }
    if( !ferr )
    {
        points[(*cnt)++] = fp;
    }

    return;
}

convexpolygon_t convexpolygon_intersect(convexpolygon_t* p, convexpolygon_t* q)
{
    int maxCnt;
    vector2_t* points;
    int cnt;
    int i;

    maxCnt = p->mPointsCnt + q->mPointsCnt;
    points = (vector2_t*)malloc(sizeof(vector2_t)*(maxCnt));
    cnt = p->mPointsCnt;

    for( i = 0; i < cnt; ++i )
    {
        points[i] = p->mPoints[i];
    }

    for( i = 0; i < q->mPointsCnt; ++i )
    {
        IntersectHalfPlane( points, &cnt, q->mPoints[i], q->mVecs[i] );
    }

    if( cnt <= 0 || q->mPointsCnt <= 0 )
    {
        free( points );
        return convexpolygon_init( 0, 0 );
    }
    else
    {
        return convexpolygon_init( points, cnt );
    }
}

convexpolyhedron_t convexpolyhedron_init
(
 vector3_t* points,
 int pointsCnt,
 edge_t* edges,
 int edgesCnt,
 face_t* faces,
 int facesCnt
)
{
    int i;

    convexpolyhedron_t retval;
    retval.mOrigPoints = malloc(sizeof(vector3_t)*pointsCnt);
    for(i = 0; i < pointsCnt; ++i) 
    {
        retval.mOrigPoints[i] = points[i];
    }
    retval.mPoints = points;
    retval.mPointsCnt = pointsCnt;
    retval.mEdges = edges;
    retval.mEdgesCnt = edgesCnt;
    retval.mFaces = faces;
    retval.mFacesCnt = facesCnt;

    return retval;
}

void dscan_convexpolyhedron_destroy
(
 convexpolyhedron_t* poly
)
{
    int i;
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        free(poly->mFaces[i].edges);
    }
    free(poly->mFaces);
    free(poly->mEdges);
    free(poly->mPoints);
    free(poly->mOrigPoints);
    free( poly );
}

void dscan_convexpolyhedron_transform
(
 convexpolyhedron_t* poly,
 matrix3_t mat
)
{
    int i;

    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        vector3_t* point = poly->mOrigPoints + i;
        float x = point->x;
        float y = point->y;
        float z = point->z;
        vector3_t* tgt = poly->mPoints + i;
        tgt->x = mat.d[0]*x + mat.d[1]*y + mat.d[2]*z;
        tgt->y = mat.d[3]*x + mat.d[4]*y + mat.d[5]*z;
        tgt->z = mat.d[6]*x + mat.d[7]*y + mat.d[8]*z;
    }
}

int dscan_convexpolyhedron_edges( convexpolyhedron_t* poly )
{
    return poly->mEdgesCnt;
}

float dscan_convexpolyhedron_volume( convexpolyhedron_t* poly )
{
    double volume = 0.0;

    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    
    vector3_t o, e0, e1, h, n;
    bool* reverse;
    int reverseCnt;

    int p, f, s;
    for( p = 0; p < poly->mPointsCnt; ++p )
    {
        vector3_t point = poly->mPoints[p];
        sx += point.x;
        sy += point.y;
        sz += point.z;
    }
    o.x = (float) (sx / poly->mPointsCnt);
    o.y = (float) (sy / poly->mPointsCnt);
    o.z = (float) (sz / poly->mPointsCnt);

    reverseCnt = 0;
    for( f = 0; f < poly->mFacesCnt; ++f )
    {
        int edgesCnt = poly->mFaces[f].edgesCnt;
        if( edgesCnt > reverseCnt )
        {
            reverseCnt = edgesCnt;
        }
    }
    reverse = malloc(sizeof(bool)*reverseCnt);

    for( f = 0; f < poly->mFacesCnt; ++f ) {
        int sides = poly->mFaces[f].edgesCnt;
        assert( sides > 2 );

        GetReverse(poly, reverse, reverseCnt, f);
        for( s = 0; s < sides - 2; ++s )
        {
            vector3_t a = GetVertexR(poly, f, s + 0, reverse[s + 0]);
            vector3_t b = GetVertexR(poly, f, s + 1, reverse[s + 1]);
            vector3_t c = GetVertexR(poly, f, s + 2, reverse[s + 2]);
            e0 = minus3(b, a);
            e1 = minus3(c, b);
            n = cross3(e0, e1);
            h = minus3(a, o);
            volume += (1.0 / 6.0) * dot3( h, n );
        }
    }

    free( reverse );
    
    return (float) volume;
}

int dscan_convexpolyhedron_visible
(
 convexpolyhedron_t* poly,
 float* edge_points
)
{
    bool* visible = malloc(sizeof(bool)*poly->mEdgesCnt);
    
    int i, j;

    int cnt = 0;

    for( i = 0; i < poly->mEdgesCnt; ++i )
    {
        visible[i] = false;
    }
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        if( IsVisible(poly, i) )
        {
            for( j = 0; j < poly->mFaces[i].edgesCnt; ++j )
            {
                visible[poly->mFaces[i].edges[j]] = true;
            }
        }
    }

    for( i = 0; i < poly->mEdgesCnt; ++i )
    {
        if( visible[i] )
        {
            edge_points[cnt++] = poly->mPoints[poly->mEdges[i].a].x;
            edge_points[cnt++] = poly->mPoints[poly->mEdges[i].a].y;
            edge_points[cnt++] = poly->mPoints[poly->mEdges[i].b].x;
            edge_points[cnt++] = poly->mPoints[poly->mEdges[i].b].y;
	}
    }

    return cnt;
}

static void fail3()
{
    fprintf( stderr, "Hard fail in convexpolyhedron_check.\n" );
    EXIT(-1);
}

static void fail3f( int face )
{
    fprintf( stderr, "Hard fail in convexpolyhedron_check, face%d.\n", face );
    EXIT(-1);
}

static void fail3fe( int face, int edge )
{
    fprintf( stderr, "Hard fail in convexpolyhedron_check, face%d, edge%d.\n",
             face, edge );
    EXIT(-1);
}

void convexpolyhedron_check
(
 convexpolyhedron_t* poly,
 bool verbose
)
{
    int i, j, k;
    float clmx;
    int clf;
    float cvmx;
    int cvf;
    int* inc;
    int* incr;
    int* incp;
    bool eifail;
    bool pifail;

    if(!poly->mPoints)
    {
        fprintf( stderr, "Points is null.\n" );
        fail3();
        return;
    }
    if( !poly->mEdges ) 
    {
        fprintf( stderr, "Edges is null.\n" );
        fail3();
        return;
    }
    if( !poly->mFaces )
    {
        fprintf( stderr, "Faces is null.\n" );
        fail3();
        return;
    }

    if( verbose ) fprintf( stderr, "Checking for degenerate faces.\n" );
    for(i = 0; i < poly->mFacesCnt; ++i) 
    {
        if( poly->mFaces[i].edgesCnt < 3 )
        {
            fail3f(i);
            return;
        }
    }

    if( verbose ) fprintf( stderr, "Checking for discontinuous faces.\n" );
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        bool* reverse = MakeReverse(poly, i);
        for( j = 0; j < poly->mFaces[i].edgesCnt; ++j )
        {
            int jp1 = ((j + 1) % poly->mFaces[i].edgesCnt);
            int end0 = reverse[j] ? poly->mEdges[poly->mFaces[i].edges[j]].a 
                : poly->mEdges[poly->mFaces[i].edges[j]].b;
            int start1 = reverse[jp1] ? poly->mEdges[poly->mFaces[i].edges[jp1]].b 
                : poly->mEdges[poly->mFaces[i].edges[jp1]].a;
            if( end0 != start1 )
                fail3fe(i, j);
        }
        free(reverse);
    }

    if( verbose ) fprintf( stderr, "Checking face colinearity:" );
    clmx = 0.0F;
    clf = -1;
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        vector3_t a = GetVertex(poly, i, 0);
        vector3_t norm = GetNorm(poly, i);

        for( j = 0; j < poly->mFaces[i].edgesCnt  - 1; ++j )
        {
            for( k = 0; k < 2; ++k )
            {
                int idx = k == 0 
                    ? poly->mEdges[poly->mFaces[i].edges[j]].a 
                    : poly->mEdges[poly->mFaces[i].edges[j]].b;
                vector3_t point = poly->mPoints[idx];
                vector3_t rpoint = minus3(point, a);
                vector3_t result = intersect_plane(norm, rpoint, norm);
                result = minus3(result, rpoint);
                if( norm3( result ) > clmx )
                {
                    clmx = norm3( result );
                    clf = i;
                }
              }
        }
    }
    if( verbose ) fprintf( stderr, " %f from Face%d\n", clmx, clf );
    /* assumes unit sphere scaling */
    if( clmx > 1e-5 )
    {
        fail3f( clf );
    }

    if( verbose ) fprintf( stderr, "Checking convexivity:" );
    cvmx = 0.0F;
    cvf = -1;
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        int index = GetVertexIndex(poly, i, 0);
        vector3_t a = GetVertex(poly, i, 0);
        vector3_t norm = GetNorm(poly, i);

        for( k = 0; k < poly->mPointsCnt; ++k )
        {
            vector3_t rel;
            float d;
            if( k == index ) continue;

            rel = minus3(poly->mPoints[k], a);
            rel = normalize3( rel );
            d = dot3( rel, norm );
            if(d > cvmx)
            {
                cvmx = d;
                cvf = i;
            }
        }
    }
    if( verbose ) fprintf( stderr, " %f from Face%d\n", cvmx, cvf );
    if( cvmx > 1e-5 )
    {
        fail3f( cvf );
    }

    if( verbose ) fprintf( stderr, "Checking edge inclusion.\n" );
    inc = malloc(sizeof(int)*poly->mEdgesCnt);
    incr = malloc(sizeof(int)*poly->mEdgesCnt);
    for( i = 0; i < poly->mEdgesCnt; ++i )
    {
        inc[i] = 0;
        incr[i] = 0;
    }
    for( i = 0; i < poly->mFacesCnt; ++i )
    {
        bool* reverse = MakeReverse(poly, i);
        for( j = 0; j < poly->mFaces[i].edgesCnt; ++j )
        {
            int idx = poly->mFaces[i].edges[j];
            if(reverse[j])
                ++incr[idx];
            else
                ++inc[idx];
        }
        free(reverse);
    }
    eifail = false;
    for( i = 0; i < poly->mEdgesCnt; ++i )
    {
        if (1 != inc[i] || 1 != incr[i]) {
            fprintf( stderr, "Edge%d has face count of %d\n",
                     i, inc[i] + incr[i] );
            eifail = true;
        }
    }
    free(inc);
    free(incr);
    if(eifail)
    {
        fail3();
    }

    if( verbose ) fprintf( stderr, "Checking point inclusion.\n" );
    incp = malloc(sizeof(int)*poly->mPointsCnt);
    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        incp[i] = 0;
    }
    for( i = 0; i < poly->mEdgesCnt; ++i )
    {
        ++incp[poly->mEdges[i].a];
        ++incp[poly->mEdges[i].b];
    }
    pifail = false;
    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        if(incp[i] <= 0)
        {
            fprintf( stderr, "Point%d has no edges\n", i );
            pifail = true;
        }
    }
    free(incp);
    if(pifail)
    {
        fail3();
    }

    if( verbose ) fprintf( stderr, "Checking Euler Formula.\n" );
    if( poly->mEdgesCnt != poly->mFacesCnt + poly->mPointsCnt - 2) {
        fprintf( stderr, "E=%d\n", poly->mEdgesCnt );
        fprintf( stderr, "F=%d\n", poly->mFacesCnt );
        fprintf( stderr, "V=%d\n", poly->mPointsCnt );
        fail3();
    }
}
