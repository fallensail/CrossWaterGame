#pragma once
#include "Graph.hpp"
#include "Point.hpp"
#include <cmath>
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

 /** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_)
      : h(h_), hx(hx_), hy(hy_) {
  }
  // More operators
  QVar& operator+=(const QVar& b) {
    h += b.h;
    hx += b.hx;
    hy += b.hy;
    return *this;
  }
  QVar& operator-=(const QVar& b) {
    h -= b.h;
    hx -= b.hx;
    hy -= b.hy;
    return *this;
  }
  QVar& operator*=(double b) {
    h *= b;
    hx *= b;
    hy *= b;
    return *this;
  }
  QVar& operator/=(double b) {
    h  /= b;
    hx /= b;
    hy /= b;
    return *this;
  }

};

QVar operator-(const QVar& a) {
  return QVar(-a.h, -a.hx, -a.hy);
}
QVar operator+(QVar a, const QVar& b) {
  return a += b;
}
QVar operator-(QVar a, const QVar& b) {
  return a -= b;
}
QVar operator*(QVar a, double b) {
  return a *= b;
}
QVar operator*(double b, QVar a) {
  return a *= b;
}
QVar operator/(QVar a, double b) {
  return a /= b;
}
/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {
  // HW3: YOUR CODE HERE
  // Write all typedefs, public member functions, private member functions,
  //   inner classes, private member data, etc.hello


 private:
  struct internal_node_value;
  struct internal_edge_value;
  struct internal_tri_value;

 public:
  /** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
  typedef unsigned size_type;
  typedef Graph<internal_node_value, internal_edge_value> GraphType;
  class Node;
  typedef Node node_type;
  class Edge;
  typedef Edge edge_type;
  class NodeIterator;
  typedef NodeIterator node_iterator;
  class EdgeIterator;
  typedef EdgeIterator edge_iterator;
  typedef N node_value_type;
  typedef E edge_value_type;
  typedef T triangle_value_type;
  /** Predeclaration of Triangle type. */
  class Triangle;
  /** Type of triangle iterators, which iterate over all adjacent mesh triangles of a triangle. */
  typedef Triangle triangle_type;
  class IncidentIterator_Triangle;  
  /** Synonym for IncidentIterator_Triangle */
  typedef IncidentIterator_Triangle incidentiterator_triangle;

  class IncidentIterator_Node;
  typedef IncidentIterator_Node incidentiterator_node;
  class IncidentEdgeIterator;
  typedef IncidentEdgeIterator incident_edge_iterator;



  Mesh() {}

  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return graph_.num_nodes();
  }

  /** Return the number of edges in the mesh. */
  size_type num_edges() const {
    return graph_.num_edges();
  }

  /** Return the number of triangles in the mesh. */
  size_type num_triangles() const {
    return tri_vec.size();
  }

  class Node:private totally_ordered<Node>{
  public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Mesh::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node(){}
    /** Return this node's position. */
    const Point& position() const {
      return mesh_->graph_.node(uid_).position();
    }

    Point& position(){
      return mesh_->graph_.node(uid_).position();
    }

    /** Return this node's index, a number in the range [0, graph_size()). */
    size_type index() const {
      return uid_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      if (x.mesh_ == mesh_ && x.uid_ == uid_)
        return true;
      return false;
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
      if (mesh_ == x.mesh_)
        return (uid_ < x.uid_);
      else
        return (mesh_ < x.mesh_);
    }

    /** Return the user-specifies value of this Node
     *  For rvalue operations
     */
    node_value_type& value(){
      return mesh_->graph_.node(uid_).value().value;
    }

    /** Return the user-specifies value of this Node
     *  For lvalue operations
     */
    const node_value_type& value() const {
      return mesh_->graph_.node(uid_).value().value;
    }

    /** Return the number of incident edges of this node
     *  Incident are edges spawned by this node
     */
    size_type degree() const{
      return mesh_->graph_.node(uid_).degree();
    }

    /* Return an iterator that points to the first incident edge of this node
     * If degree()==0, the returned iterator value shall not be dereferenced
     */
    incident_edge_iterator edge_begin() const{
      return IncidentEdgeIterator(mesh_, mesh_->graph_.node(uid_).edge_begin());
    }
    /*  Returns an iterator referring to the past-the-end incident edge
       *  of this node
     *The past-the-end incident edge is the theoretical incident edge that 
       *  would follow the last edge. It does not point to any element, and 
       *  thus shall not be dereferenced.
     *If degree()==0, this function returns the same as edge_begin().
     */
    incident_edge_iterator edge_end() const{
      return IncidentEdgeIterator(mesh_, mesh_->graph_.node(uid_).edge_end());
    }

    /* Get the iterator that points to the first adjacent triangle of a node
     * @param n The node in the center
     * @return the iterator that points to the first adjacent triangle of @a n
     * @post result.n_ == n.uid_
     * @post result.t_ result.t2_ are the adjacent triangles of edge(n_, result.last_)
     * @post result.last_ == result.first_
     */
    IncidentIterator_Node triangle_begin(){
      Edge start = *(edge_begin());
      size_type othernode;
      if (start.node1().index() == index())
        othernode = start.node2().index();
      else
        othernode = start.node1().index();
      return IncidentIterator_Node(mesh_, index(), mesh_->graph_.edge(start.index()).value().triangle1, 
        mesh_->graph_.edge(start.index()).value().triangle2, othernode, othernode);
    }

    /* Get the iterator that points to an invalid triangle of a node
     * @param n The node in the center
     * @return the iterator that points to an invalid triangle of @a n
     * @post result.n_ == n.uid_
     * @post result.t_ == -1
     */
    IncidentIterator_Node triangle_end(){
      return IncidentIterator_Node(mesh_, index(), -1, 
        -1, -1, -1);
    }

  private:
    friend class Mesh;
    Mesh* mesh_; //pointer to the associated Mesh
    size_type uid_; //uid of the node; the same uid as in graph
    // a private constructor for graph to construct a Node instance
    Node(const Mesh* mesh, size_type uid)
      : mesh_(const_cast<Mesh*>(mesh)), uid_(uid){
      }
  };

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value_in The value associated with this node
   * @post new graph_.size() == old graph_.size() + 1
   * @post result_node.index() == old size()
   * Complexity: O(1)
   */
  Node add_node(const Point& position, const N& value_in = N()) {
    auto n = graph_.add_node(position, internal_node_value(value_in));
    return Node(this, n.index());
  }

  Node node(size_type i){
    return Node(this, i);
  }

  class Edge:private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return mesh_->node(mesh_->graph_.edge(uid_).node1().index());
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return mesh_->node(mesh_->graph_.edge(uid_).node2().index());
    }

    double length() const{
      return mesh_->graph_.edge(uid_).length();
    }

    Triangle triangle1() {
      return Triangle(mesh_, mesh_->graph_.edge(uid_).value().triangle1);
    }

    Triangle triangle2() {
      return Triangle(mesh_, mesh_->graph_.edge(uid_).value().triangle2);
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      if (mesh_ == x.mesh_ && uid_ == x.uid_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
      if (mesh_ == x.mesh_)
        return (uid_ < x.uid_);
      else
        return (mesh_ < x.mesh_);
    }
    
    //Return the user-specifies value of this edge
    edge_value_type& value(){
      return mesh_->graph_.edge(uid_).value().value;
    }

    //Return the user-specifies value of this edge
    const edge_value_type& value() const{
      return mesh_->graph_.edge(uid_).value().value;
    }

    size_type index() const{
      return uid_;
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;  //pointer to the associated graph
    size_type uid_; //the unique id of the edge
    
    //private constructor for graph to construct edge instance
    Edge(const Mesh* mesh, size_type uid)
      : mesh_(const_cast<Mesh*>(mesh)), uid_(uid){
      }
  };

  Edge edge(size_type i){
    return Edge(this, i);
  }

  class Triangle
  {
  public:
    /**construc an invalid triangle*/
    Triangle():uid_(-1){
    }
    /**Access the i node of the triangle
     * @pre 0<=i<3
     * @post result.index() == mesh_->tri_vec[uid_].nodes[i]
     * Complexity = O(1)
     */
    Node node(size_type i) {
      size_type node_uid = mesh_->tri_vec[uid_].nodes[i];
      return mesh_->node(node_uid);
    }
    /**Access the first node of the triangle
     * @pre 0<=i<3
     * @post result.index() == mesh_->tri_vec[uid_].edges[i]
     * Complexity = O(1)
     */
    Edge edge(size_type i) {
      size_type edge_uid = mesh_->tri_vec[uid_].edges[i];
      return mesh_->edge(edge_uid);
    }

    /**return the area of this triangle
     *Comlexity = O(1)
     */
    double area() {
      mesh_->tri_vec[uid_].area = norm(cross(node(0).position()-node(1).position(), 
        node(0).position()-node(2).position()))/2;
      return mesh_->tri_vec[uid_].area;
    }
    /**return the index of the triangle
     *Complexity = O(1)
     */
    size_type index(){
      return uid_;
    }

    /** Return the value associated with this Node
     *  For rvalue operations
     */
    triangle_value_type& value(){
      return mesh_->tri_vec[uid_].value;
    }

    /** Return the value associated with this Node
     *  For lvalue operations
     */
    const triangle_value_type& value() const {
      return mesh_->tri_vec[uid_].value;
    }



    /**return the unit normal vector of the edge
        *that is opposite to node i
      */
    Point normal(size_type i){
      return mesh_->tri_vec[uid_].n[i];
    }

    /* Return an iterator points to the first adjacent triangle of this triangle.
     * @post result->incident_i ==0
     * @post result->mesh_ == this->mesh_
     * @post result->triangle_uid == this->uid_
     *
     *Complexity = O(1)
         */
    incidentiterator_triangle triangle_begin() const {
      return IncidentIterator_Triangle(mesh_, uid_, 0);
    }
  
        /* Return an iterator points to the last adjacent triangle of this triangle.
     * @post result->incident_i == this->mesh_->edges.size()
     * @post result->mesh_ == this->mesh_
     * @post result->triangle_uid == this->uid_
     *
     *Complexity = O(1)
         */
    incidentiterator_triangle triangle_end() const {
      return IncidentIterator_Triangle(mesh_, uid_, 3);
    }

  private:
    friend class Mesh;
    Triangle(Mesh* mesh, size_type uid):
      uid_(uid), mesh_(mesh) {}
    size_type uid_;
    Mesh* mesh_;
  };



  /** Add an triangle to the graph
   *@param[in] a, b, c are all valid nodes in the trigraph
   *@return The triangle that has been added
   *@post new size() = old size() +1
   * Complexity: O(d) d is the degree of a node
   */
  Triangle add_triangle(const Node& a, const Node& b, const Node& c){
    internal_triangle_value new_triangle;
    //add edge in graph
    auto edge0 = graph_.add_edge(graph_.node(a.index()), graph_.node(b.index()));    
    new_triangle.nodes[0] = c.index();
    new_triangle.edges[0] = edge0.index();
    
    auto edge1 = graph_.add_edge(graph_.node(b.index()), graph_.node(c.index()));    
    new_triangle.nodes[1] = a.index();
    new_triangle.edges[1] = edge1.index();
    
    auto edge2 = graph_.add_edge(graph_.node(a.index()), graph_.node(c.index()));    
    new_triangle.nodes[2] = b.index();
    new_triangle.edges[2] = edge2.index();    
    double xa = a.position().x;
    double ya = a.position().y;
    double xb = b.position().x;
    double yb = b.position().y;
    double xc = c.position().x;
    double yc = c.position().y;
    new_triangle.area = (xa*yb+xb*yc+xc*ya-xb*ya-xc*yb-xa*yc)/2;
    if (new_triangle.area < 0)
      new_triangle.area = -new_triangle.area;

    Point p = a.position() - b.position(); 
    double nx = p.y;
    double ny = 0-p.x;
    // check if it is in the right direction
    Point checkp = c.position() - a.position();
    // check direction, if wrong, flip it. & adjust length to be equal to |e|
    if(nx*checkp.x + ny*checkp.y>0) {
      nx = -nx;
      ny = -ny;
    }  
    new_triangle.n.push_back(Point(nx, ny, 0));
    p = b.position() - c.position(); 
    nx = p.y;
    ny = 0-p.x;
    
    // check if it is in the right direction
    checkp = a.position() - b.position();
    // check direction, if wrong, flip it. & adjust length to be equal to |e|
    if(nx*checkp.x + ny*checkp.y>0) {
      nx = -nx;
      ny = -ny;
    }  
    
    new_triangle.n.push_back(Point(nx, ny, 0));     
    p = a.position() - c.position(); 
    nx = p.y;
    ny = 0-p.x;
    // check if it is in the right direction
    checkp = b.position() - a.position();
    // check direction, if wrong, flip it. & adjust length to be equal to |e|
    if(nx*checkp.x + ny*checkp.y>0) {
      nx = -nx;
      ny = -ny;
    }  
    new_triangle.n.push_back(Point(nx, ny, 0));
      
    tri_vec.push_back(new_triangle);
    
    if(edge0.value().triangle1 == unsigned(-1))
      edge0.value().triangle1 = tri_vec.size()-1;
    else
      edge0.value().triangle2 = tri_vec.size()-1;
    if(edge1.value().triangle1 == unsigned(-1))
      edge1.value().triangle1 = tri_vec.size()-1;
    else
      edge1.value().triangle2 = tri_vec.size()-1;
    if(edge2.value().triangle1 == unsigned(-1))
      edge2.value().triangle1 = tri_vec.size()-1;
    else
      edge2.value().triangle2 = tri_vec.size()-1;
      return Triangle(this,tri_vec.size()-1);
  }

  /** Determine if this Triangle belongs to this Graph
    * @return True if @a t is currently a Triangle of this Graph
    *
    * Complexity: O(1).
    */
  bool has_triangle(const Triangle& t){
    //compare the trigraph_ and uid_
    return ((t.graph_ != this->graph_) && (t.uid_ < tri_vec.size()));   
  }

    /** Return the triangle with index @a i.
     * @pre 0 <= @a i < size()
     * @post result.index() == i
     * Complexity: O(1).
     */
  Triangle triangle(size_type i){
    assert(i<tri_vec.size());
    return Triangle(this,i);
  }
  
  
  // /** Return the first triangle adjacent to the input edge.
  //    * @param[in] e, an edge inside the mesh
  //    * @pre @a e.index() < graph_->num_edges()
  //    * @post result.uid_ == graph_->get_edge_value(edge_uid).F_triangle1.
  //    * Complexity: O(1).
  //    */
  // Triangle triangle1(const Edge& e){
  //   return Triangle(this,graph_.edges[e.index()].triangle1);
  // }

  // /** Return the second triangle adjacent to the input edge.
  //    * @param[in] e, an edge inside the mesh
  //    * @pre @a e.index() < graph_->num_edges()
  //    * @post result.uid_ == graph_->get_edge_value(edge_uid).F_triangle2.
  //    * Complexity: O(1).
  //    */
  // Triangle triangle2(const Edge& e){
  //   return Triangle(this,graph_.edges[e.index()].triangle2);
  // }

  class NodeIterator:private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    //return the Node pointed by the interator
    Node operator*() const{
      return mesh_->node((*it_).index());
    }

    NodeIterator& operator++(){
      ++it_;
      return *this;
    }

    bool operator==(const NodeIterator& x) const{
      return (it_ == x.it_);
    }

  private:
    friend class Mesh;
    Mesh* mesh_;
    typename GraphType::node_iterator it_;
    //private constructor for graph to construct NodeIterator instance
    NodeIterator(const Mesh* mesh, typename GraphType::node_iterator it)
      : mesh_(const_cast<Mesh*>(mesh)), it_(it){
      }

  };


  /* Return an iterator points to the first node in graph_.
     * Complexity: O(1).
     */
  node_iterator  node_begin(){
    return NodeIterator(this, graph_.node_begin());
  }

    /* Return an iterator points to the last node in graph_.
     * Complexity: O(1).
     */
  node_iterator node_end(){
    return NodeIterator(this, graph_.node_end());
  }

  class EdgeIterator:private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    //return the Edge pointed by the interator
    Edge operator*() const{
      return mesh_->edge((*it_).index());
    }

    EdgeIterator& operator++(){
      ++it_;
      return *this;
    }

    bool operator==(const EdgeIterator& x) const{
      return (it_ == x.it_);
    }

  private:
    friend class Mesh;
    Mesh* mesh_;
    typename GraphType::edge_iterator it_;
    //private constructor for graph to construct EdgeIterator instance
    EdgeIterator(const Mesh* mesh, typename GraphType::edge_iterator it)
      : mesh_(const_cast<Mesh*>(mesh)), it_(it){
      }
  };


    
    /* Return an iterator points to the first edge in graph_.
     * Complexity: O(1).
     */
  edge_iterator edge_begin(){
    return EdgeIterator(this, graph_.edge_begin());
  }

    /* Return an iterator points to the last edge in graph_.
     * Complexity: O(1).
     */
  edge_iterator edge_end(){
    return EdgeIterator(this, graph_.edge_end());
  }

  class IncidentEdgeIterator:private totally_ordered<IncidentEdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentEdgeIterator. */
    IncidentEdgeIterator() {
    }

    //return the Edge pointed by the interator
    Edge operator*() const{
      return mesh_->edge((*it_).index());
    }

    IncidentEdgeIterator& operator++(){
      ++it_;
      return *this;
    }

    bool operator==(const IncidentEdgeIterator& x) const{
      return (it_ == x.it_);
    }

  private:
    friend class Mesh;
    Mesh* mesh_;
    typename GraphType::incident_iterator it_;
    //private constructor for graph to construct EdgeIterator instance
    IncidentEdgeIterator(const Mesh* mesh, typename GraphType::incident_iterator it)
      : mesh_(const_cast<Mesh*>(mesh)), it_(it){
      }
  };


  //Iterator that iterates through the adjacent triangle of a node
  class IncidentIterator_Node:private totally_ordered<IncidentIterator_Node>
  {
  public:
    // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Triangle value_type;
      /** Type of pointers to elements. */
      typedef Triangle* pointer;
      /** Type of references to elements. */
      typedef Triangle& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;
    

    /** Construct an invalid IncidentIterator_Node*/
    IncidentIterator_Node(){
    }

    /* Return the Triangle pointed by the IncidentIterator_Node
     * @post return.uid_ = t_
     */
    Triangle operator*() const{
      return mesh_->triangle(t_);
    }

    /* return the IncidentIterator_Node that points to the
     * next incident triangle
     * @pre t_ old != -1
     * @post If t_ old has two adjacent triangles, init_ != last_ 
       * t_ new points to the adjacent triangle which does not 
       * contain node last_.
       * If t_ old has two adjacent triangles, init_ == last_, t_ new = -1
       * If t_ old has one adjacent triangles, t2_!=-1, t_ new = t2_
       * If t_ old has one adjacent triangles, t2_ == -1, t_ new = -1
     */
    IncidentIterator_Node& operator++(){
      //get an incident edge of the node
      //record the last_ = first_ = incident_edge.node().index()
      //t_ = incident_edge.value().triangle1 t2_ = incident_edge.value().triangle2
      //get the new t_ and update last_
      //if t_ == t2_ set t_ to invalid to mark the end
      //if meets the border. set last_ = init_ and explore t2_
      //when t2_ meets the end, set t_ to invilid

      if (t_ == t2_){
        t_ = -1;
      }
      else{
        Triangle old = mesh_->triangle(t_);
        size_type i = 0;
        while(old.node(i).index() != last_)
          ++i;
        if (old.mesh_->graph_.edge(old.edge(i).index()).value().triangle1 == unsigned(-1) || 
          old.mesh_->graph_.edge(old.edge(i).index()).value().triangle2 == unsigned(-1)){
          if (t2_ == unsigned(-1))
            t_ = -1;
          else{
            t_ = t2_;
            t2_ = -1;
            last_ = init_;
          }
        }
        else{
          if (old.edge(i).node1().index() == n_)
            last_ = old.edge(i).node2().index();
          else
            last_ = old.edge(i).node1().index();
          if (old.mesh_->graph_.edge(old.edge(i).index()).value().triangle1 == t_)
            t_ = old.mesh_->graph_.edge(old.edge(i).index()).value().triangle2;
          else
            t_ = old.mesh_->graph_.edge(old.edge(i).index()).value().triangle1;
        }
      }
      return *this;
      }

      //True if n_ and t_ are the same
      bool operator==(const IncidentIterator_Node& x) const{
        if (mesh_ == x.mesh_ && t_ == x.t_ && n_ == x.n_)
          return true;
        return false;
      }
  private:
    friend class Mesh;
    Mesh* mesh_;
    size_type n_; //The uid of the node in the center
    size_type t_; //The uid of the triangle pointed to
    size_type t2_; //the uid of the triangle to be explored when t_ reaches boundary
    size_type last_; //the uid of the other node of the edge just discovered
    size_type init_; //the uid of the other node of the first edge

    /** Construct a valid IncidentIterator_Node*/
    IncidentIterator_Node(const Mesh* m, size_type n, size_type t, size_type t2, size_type last, size_type init): mesh_(const_cast<Mesh*>(m)),n_(n),t_(t),t2_(t2),last_(last),init_(init){
    }
  };



  /* Get the iterator that points to the first adjacent triangle of a node
   * @param n The node in the center
   * @return the iterator that points to the first adjacent triangle of @a n
   * @post result.n_ == n.uid_
   * @post result.t_ result.t2_ are the adjacent triangles of edge(n_, result.last_)
   * @post result.last_ == result.first_
   */
  IncidentIterator_Node triangle_begin(Node n){
    Edge start = *(n.edge_begin());
    size_type othernode;
    if (start.node1().index() == n.index())
      othernode = start.node2().index();
    else
      othernode = start.node1().index();

    
    return IncidentIterator_Node(this, n.index(), graph_.edge(start.index()).value().triangle1, 
      graph_.edge(start.index()).value().triangle2, othernode, othernode);

  }

  /* Get the iterator that points to an invalid triangle of a node
   * @param n The node in the center
   * @return the iterator that points to an invalid triangle of @a n
   * @post result.n_ == n.uid_
   * @post result.t_ == -1
   */
  IncidentIterator_Node triangle_end(Node n){
    return IncidentIterator_Node(this, n.index(), -1, 
      -1, -1, -1);
  }


  //Iterators:
  //all triangles
  //adjacent triangles of a node
  //adjacent triangles of a triangle

  //Iterator that iterates through the adjacent triangle of a triangle
  /** @class Mesh::IncidentIterator_Triangle
     * @brief Interator Class for triangles incident to a triangle. A forward iterator.
     * @RI graph_pointer != nullptr && incident_i <= 2
     */
  class IncidentIterator_Triangle:private totally_ordered<IncidentIterator_Triangle>
  {
  public:
    // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Triangle value_type;
      /** Type of pointers to elements. */
      typedef Triangle* pointer;
      /** Type of references to elements. */
      typedef Triangle& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;
    

    /** Construct an invalid IncidentIterator_Triangle*/
      IncidentIterator_Triangle(){
      }

      /**return the Triangle pointed by the IncidentIterator_Triangle
       *@pre this != nullptr
       *@post result is a triange adjacent to Triangle(this->mesh_,this->triangle_uid)
       *@post result is one of the two adjacent triangles to the ith edge of Triangle(this->mesh_,this->triangle_uid)  
       *@post result.uid_ != this->triangle_uid
       *
       *Complexity = O(1)
       */
      Triangle operator*() {
        size_type e_uid = mesh_->tri_vec[triangle_uid].edges[incident_i];

       
        //Edge e = mesh_->graph_.edge(e_uid);
        if(mesh_->graph_.edge(e_uid).value().triangle1 != triangle_uid)
          return Triangle(mesh_,mesh_->graph_.edge(e_uid).value().triangle1);
        else 
          return Triangle(mesh_,mesh_->graph_.edge(e_uid).value().triangle2);
        
       
      }

      /*return the IncidentIterator_Triangle that points to the next incident edge
       *@pre this->incident_i < 2
       *@post result.incident_i > this->incident_i
       *@post result.incident_i <=2
       *
       *Complexity = O(1)
       */
      IncidentIterator_Triangle& operator++(){
        incident_i++;
        return *this;
      }

        /*return true if the two triangles are same and false if not
         *@param[in] x, reference to a triangle 
       *@pre @a x.uid_ >=0
       *@pre triangle_uid >=0
       *@post result==true if @a x.uid_ == this->triangle_uid  &&  @a x.mesh_ == this->mesh_
       *@post result==false if @a x.uid_ != this->triangle_uid ||  @a x.mesh_ != this->mesh_
       *
       *Complexity = O(1)
       */
      bool operator==(const IncidentIterator_Triangle& x) const{
        return((x.triangle_uid == this->triangle_uid)  &&  (x.mesh_ == this->mesh_) && (x.incident_i == this->incident_i));
      }
      
      size_type index(){
        return incident_i;
      }
  private:
    friend class Mesh;
    Mesh* mesh_;  //poiter pointing to this Mesh
    size_type triangle_uid;  //the uid of this triangle
    /** incident_i represents the ith edge of this triangle*/
    size_type incident_i;

        /** Construct a valid IncidentIterator_Triangle*/
        IncidentIterator_Triangle(const Mesh* m, size_type uid, size_type i): mesh_(const_cast<Mesh*>(m)),triangle_uid(uid),incident_i(i){
      }
  };

  class TriangleIterator: private totally_ordered<TriangleIterator>{
  public:

    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    TriangleIterator(){}

    Triangle operator*() const{
      return mesh_->triangle(idx_);
    }

    TriangleIterator& operator++(){
      ++idx_;
      return *this;
    }

    bool operator==(const TriangleIterator& x) const{
      if (mesh_ == x.mesh_ && idx_ == x.idx_)
        return true;
      return false;
    }

  private:
    friend class Mesh;

    TriangleIterator(const Mesh* mesh, size_type idx)
      : mesh_(const_cast<Mesh*>(mesh)), idx_(idx){
      }
    Mesh* mesh_;
    size_type idx_;
  };

  TriangleIterator triangle_begin() const{
    return TriangleIterator(this, 0);
  }

  TriangleIterator triangle_end() const{
    return TriangleIterator(this, tri_vec.size());
  }


private:

  struct internal_triangle_value
  {
    triangle_value_type value;
    std::vector<size_type> nodes; //a vector storing the uids of the three nodes of this triangle
    std::vector<size_type> edges; //a vector storing the uids of the three edges of this triangle
    double area;  //the area of this triangle
    std::vector<Point> n; //The unit normal vector of 3 edges
    internal_triangle_value():value(triangle_value_type()), nodes(3,0),
        edges(3,0), area(-1){}
  };

  struct internal_node_value
  {
    node_value_type value;
    internal_node_value() {}
    internal_node_value(const node_value_type& v):value(v){}
  };

  struct internal_edge_value
  {
    edge_value_type value;
    size_type triangle1;
    size_type triangle2;
    internal_edge_value():triangle1(-1), triangle2(-1){}
  };

  GraphType graph_; //a Graph Class object storing all the nodes and edges
  std::vector<internal_triangle_value> tri_vec; // a vector storing all the triangles' data
};
