/***************************************************************************
* abvh.cpp   (aggregate object plugin)                                     *
*                                                                          *
* The abvh (Automatic Bounding Volume Hierarchy) creates a bounding        *
* volume hierarchy using branch-and-bound to find an optimal hierarchy     *
* with respect to a probabilistic objective function and the given         *
* ordering of the objects.  That is, the hierarchy depends upon the order  *
* in which the objects are inserted.                                       *
*                                                                          *
* History:                                                                 *
*   10/09/2005  Ported from a previous ray tracer.                         *
*                                                                          *
***************************************************************************/
#include "ray_tracer.h"
#include "util.h"
#include "params.h"


// The node struct forms all of the nodes in the bounding volume hierarchy,
// which can be a tree with arbitrary branching (i.e. internal nodes can
// have any number of children).
struct node {
    node();
   ~node() {}
    double AdjustedInternalCost() const;
    void   AddChild( node * );
    bool   Leaf() const { return child == NULL; }
    const Object *object;
    AABB    bbox;    // The bounding box containing all the children.
    double  EC;      // External cost of this object.
    double  SA;      // Surface area -- can be zero for primitve.
    double  AIC;     // Adjusted internal cost.
    double  SEC_;    // Sum of external costs of children.
    double  SAIC;    // Sum of adjusted internal costs of children.
    node   *sibling; // The next child of the node that contains this one. 
    node   *parent;  // The node one level up, or NULL if this is the root.
    node   *child;   // An object or volume nested inside this one.
    };

struct abvh : public Aggregate { 
    abvh() { abvhroot = NULL; }
   ~abvh() {}
    virtual bool Intersect( const Ray &ray, HitInfo & ) const;
    virtual bool Inside( const Vec3 & ) const;
    virtual Interval GetSlab( const Vec3 & ) const;
    virtual string MyName() const { return "abvh"; }
    virtual void Close();
    virtual double Cost() const;
    static bool Branch_and_Bound( node*, node*, node*&, double& ); 
    void Insert( const Object *, double relative_cost = 1.0 );
    node *abvhroot;
    };

// When the object is closed, add each of the child objects to the hierarchy.
// Waiting until we have all the objects allows us to (otionally) radomize the
// order of insertion, which can greatly affect the quality of the resulting
// bounding volume hierarchy.
void abvh::Close()
    {
    // Should "randomize" here...
    for( unsigned i = 0; i < NumChildren(); i++ )
        {
        const Object *obj = GetChild(i);
        // Take into account the estimated cost of intersecting a ray with
        // this object as we insert it into the existing bvh.
        Insert( obj, obj->Cost() );
        }
    }


// The node_iterator struct is the mechanism for walking through a bounding
// volume hierarchy, either accessing all the volumes (in-order traversal),
// or skipping sub-trees.
struct node_iterator {
    node_iterator( const node *r ) { curr = r; }
   ~node_iterator() {}
    inline bool IsLeaf() const { return curr->child == NULL; }
    inline const node *SkipChildren(); // Next volume, ignoring children.
    inline const node *Next(); // Down if possible, then to sibling, or up.
    inline const node *Curr() const { return curr; }
    private: const node *curr;
    };

inline const node *node_iterator::SkipChildren()
    {
    // Advance the "current" pointer (curr) to the next unvisited node
    // at the same level or higher.  Since the parent has already been
    // visited, when we move "up" we must also move to the next "sibling".
    while( curr != NULL )
        {
        if( curr->sibling != NULL ) { curr = curr->sibling; break; }
        curr = curr->parent;
        }
    return curr;
    }

inline const node *node_iterator::Next()
    {
    // Find the next child in the hierarchy by first going "down" to
    // a child node, if possible.  When there is no child node, try
    // moving to a sibling.  If that fails, move up to the parent first
    // and then to a sibling (i.e. try to find an uncle).
    if( curr->child == NULL ) SkipChildren();
    else curr = curr->child;
    return curr;
    }

bool abvh::Intersect( const Ray &ray, HitInfo &hitinfo ) const
    {
    // Walk the bounding volume hierarchy intersecting, descending down a
    // branch if and only if the ray intersects the volume.  When a leaf is
    // reached, intersect the ray with the object found there.
    bool found_a_hit = false;
    double closest_hit = Infinity;
    node_iterator iter( abvhroot );
    while( iter.Curr() != NULL )
        {
        if( iter.IsLeaf() )
            {
            if( iter.Curr()->object->Intersect( ray, hitinfo ) )
                {
                found_a_hit = true;
                closest_hit = hitinfo.distance;
                }
            iter.Next();
            }
        else if( Hit( ray, iter.Curr()->bbox, closest_hit ) ) iter.Next();
        else iter.SkipChildren();
        }
    return found_a_hit;
    }

// Test to see if the point P is "inside" the object.  For an aggregate object
// this is done by asking each of the child objects.
bool abvh::Inside( const Vec3 &P ) const
    {
    node_iterator iter( abvhroot );
    while( iter.Curr() != NULL )
        {
        if( iter.IsLeaf() && iter.Curr()->object->Inside( P ) ) return true;
        iter.Next();
        }
    return false;
    }

// Return an interval that bounds the object in the given direction.  For an
// aggregate this is done by consulting all the child objects. 
Interval abvh::GetSlab( const Vec3 &v ) const
    {
    Interval I;
    node_iterator iter( abvhroot );
    while( iter.Curr() != NULL )
        {
        // For each leaf node encountered, expand the interval to
        // include the interval of the child object.
        if( iter.IsLeaf() ) I << iter.Curr()->object->GetSlab( v ); 
        iter.Next();
        }
    return I;
    }

// This function adds an object to an existing bounding volume hierarchy.
// It first constructs a bounding box for the new object, then finds the
// optimal place to insert it into the hierachy using branch-and-bound, and
// finally inserts it at the determined location, updating all the costs
// associated with the nodes in the hierarchy as a side-effect.
void abvh::Insert( const Object *obj, double relative_cost )
    {
    AABB box( GetBox( *obj ) );
    double bound = 0.0;
    node *n   = new node;
    n->bbox   = box;
    n->object = obj;
    n->SA     = SurfaceArea( box );
    n->EC     = 1.0;
    n->SEC_   = relative_cost;
    n->AIC    = relative_cost * n->SA;
    n->SAIC   = 0.0; // There are no child volumes yet.

    // If this is the first node being added to the hierarchy,
    // it becomes the root.
    if( abvhroot == NULL ) { abvhroot = n; return; }

    // Look for the optimal place to add the new node.
    // The branch-and-bound procedure will figure out where
    // to put it so as to cause the smallest increase in the
    // estimated cost of the tree.
    bound = Infinity;
    node *insert_here = NULL;
    Branch_and_Bound( abvhroot, n, insert_here, bound );
    // Now actually insert the node in the optimal place.
    insert_here->AddChild( n );
    }

// Estimate the average cost of intersecting a ray with this object.
double abvh::Cost() const
    {
    if( abvhroot == NULL ) return 1.0;
    return abvhroot->EC + abvhroot->SA * abvhroot->AIC;
    }

// Initialize all the fields of the new abvh object.
node::node()
    {
    EC      = 0.0;
    SA      = 0.0;
    AIC     = 0.0;
    SEC_    = 0.0;
    SAIC    = 0.0;
    sibling = NULL;
    parent  = NULL;
    child   = NULL;
    object  = NULL;
    bbox    = AABB::Null();
    }

 // The "Adjusted Internal Cost" of a bounding volume is the product of
 // its surface area and its internal cost.  By definition, this is zero
 // for primitive objects.
double node::AdjustedInternalCost() const
    {
    // Both "Internal Cost" and "Adjusted Internal Cost" of a
    // leaf are zero by definition.
    if( child == NULL ) return 0.0;
    double sum_EC  = 0;
    double sum_AIC = 0;
    for( const node *c = child; c != NULL; c = c->sibling )
        {
        sum_EC  += c->EC;
        sum_AIC += c->AdjustedInternalCost();
        }
    return SA * sum_EC + sum_AIC;
    }

// This method determines the optimal location to insert a new node
// into the hierarchy; that is, the place that will create the smallest
// increase in the expected cost of intersecting a random ray with the
// new hierarchy.  This is accomplished with branch-and-bound, which finds
// the same solution as a brute-force search (i.e. attempting to insert
// the new object at *each* possible location in the hierarchy), but
// does so efficiently by pruning large portions of the tree.
bool abvh::Branch_and_Bound( 
    node  *the_root,     // Root of the tree to be modified.
    node  *the_node,     // The node to add to the hierarchy.
    node  *&best_parent, // The best node to add it to.
    double &bound )      // Smallest bound achieved thus far.

    {
    double a_delta;  // AIC increase due to Root bbox expansion. 
    double r_delta;  // AIC increase due to new child of Root.   
    double c_bound;  // bound - a_delta, or less.               

    // Compute the increase in the bounding box of the root due 
    // to adding the new object.  This is used in computing     
    // the new area cost no matter where it ends up.            
    AABB root_box( the_root->bbox );
    bool expanded = root_box.Expand( the_node->bbox );

    // Compute the common sub-expression of the new area cost  
    // and the increment over this that would occur if we     
    // made the new object a child of the_root.                    
    if( expanded )
        {
        double new_area = SurfaceArea( root_box );
        a_delta  = ( new_area - the_root->SA ) * the_root->SEC_;
        c_bound  = bound - a_delta;
        if( c_bound <= 0.0 ) return false; // Early cutoff. 
        r_delta  = new_area * the_node->EC + the_node->AIC;
        }
    else
        {
        a_delta = 0.0;
        c_bound = bound;
        r_delta = the_root->SA * the_node->EC + the_node->AIC;
        }

    // See if adding the new node directly to the root of this tree
    // achieves a better bound than has thus far been obtained (via
    // previous searches).  If so, update all the parameters to reflect
    // this better choice.
    bool improved = false;
    if( r_delta < c_bound )
        {
        bound       = a_delta + r_delta;  // This is < *bound. 
        best_parent = the_root;
        c_bound     = r_delta;
        improved    = true;
        }

    // Compute the smallest increment in total area cost (over 
    // "common") which would result from adding the new node 
    // to one of the children of the_root.  If any of them obtains
    // a better bound than achieved previously, then the "best_parent"
    // and "bound" parameters will be updated as a side effect.
    for( node *c = the_root->child; c != NULL; c = c->sibling )
        {
        if( Branch_and_Bound( c, the_node, best_parent, c_bound ) )
            {
            bound   = c_bound + a_delta;
            improved = true;
            }
        }
    return improved;
    }

// This function inserts the given child node into the bounding volume
// hierarchy and adjusts all the associated costs that are stored in
// the tree.
void node::AddChild( node *new_child )
    {
    double old_AIC;
    double increment;
    AABB new_volume( new_child->bbox );

    if( child == NULL ) 
        {
        // The current node is a leaf, so we must convert it into
        // an internal node.  To do this, copy its current contents into
        // a new node, which will become a child of this one.
        child         = new node;
        child->bbox   = bbox;
        child->object = object;
        child->parent = this;
 
        // Fill in all the fields of the current node, which has just changed
        // into an internal node with a single child.
        EC   = 1;
        SA   = child->SA;
        SEC_ = child->EC;  // There is only one child.
        SAIC = child->AIC; // There is only one child.
        AIC  = SA * SEC_ + SAIC;
        }

    // Splice the new child into the linked list of children.
    new_child->sibling = child;
    new_child->parent  = this;
    child = new_child;
    
    // Update the summed external cost and the summed AIC as a result
    // of adding the new child node.  These do not depend on surface
    // area, so we needn't consider any expansion of the bounding volume.
    SEC_ += new_child->EC;
    SAIC += new_child->AIC;

    // Now take bounding volume expansion into account due to the new child.
    bool expanded = false;
    if( bbox.Expand( new_volume ) )
        {
        expanded = true;
        SA = SurfaceArea( bbox );
        }

    // Compute new area cost & how much it increased due to new child.
    old_AIC = AIC;
    AIC = SA * SEC_ + SAIC;
    increment = AIC - old_AIC;

    // Propagate the information upward, in two phases.  The first phase
    // deals with volumes that get expanded.  Once a volume is reached that does
    // not expand, there will be no more expansion all the way to the root.
    node *n = this;
    while( expanded )
        {
        expanded = false;
        node *prev = n;
        n = n->parent;
        if( n != NULL && n->bbox.Expand( new_volume ) )
            {
            expanded  = true;
            old_AIC   = n->AIC;
            n->SA     = SurfaceArea( n->bbox );
            n->SAIC  += increment;
            n->AIC    = n->SA * n->SEC_ + n->SAIC;
            increment = n->AIC - old_AIC;
            }
        }

    // From this point up to the root there will be no more expansion.
    // However, we must still propagate information upward.
    while( n != NULL )
        {
        n->SAIC += increment;
        n->AIC  += increment;
        n = n->parent;
        }
    }
    
