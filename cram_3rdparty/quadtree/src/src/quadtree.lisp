#|
  This file is a part of quadtree project.
  Copyright (c) 2015 Masayuki Takagi (kamonama@gmail.com)
|#

(in-package :cl-user)
(defpackage quadtree
  (:use :cl :cl-transforms)
  (:export #:make-quadtree-internal
           #:make-quadtree
           #:matrix->quadtree
           #:boundary
           #:insert
           #:query
           #:clear-quadtree
           #:get-in-quadtree
           #:point-intersect-p
           #:merge-quadtrees
           #:object-equal
           #:same-coords-as
           #:normalize-quadtree
           #:quadtree-map))
(in-package :quadtree)

(defvar *max-depth* 50)
(defvar *middlepoint-offset* 0.00001) ;; [mm/10]

(defstruct (quadtree (:constructor %make-quadtree)) ;; % prefix, because one should use make-quadtree instead
  object
  nw ne sw se
  boundary
  depth
  max-depth)

(defun make-quadtree (x0 y0 x1 y1 &key (max-depth *max-depth*))
  (%make-quadtree-internal x0 y0 x1 y1 0 max-depth))

;; (defmethod initialize-instance :after ((grid occupancy-grid) &key)
;;   (with-slots (width height resolution grid) grid
;;     (setf grid (make-array (list (round (/ height resolution))
;;                                  (round (/ width resolution)))
;;                            :element-type 'fixnum
;;                            :initial-element 0))))

(defun %make-quadtree-internal (x0 y0 x1 y1 depth max-depth)
  (declare (float x0 y0 x1 y1))
  (unless (<= 0 depth max-depth)
    (error "The value ~S is an linvalid value as depth." depth))
  (%make-quadtree :boundary (list x0 y0 x1 y1)
                  :depth depth
                  :max-depth max-depth))

(defun quadtree-x0 (quadtree-boundary)
  (when quadtree-boundary
    (first quadtree-boundary)))

(defun quadtree-y0 (quadtree-boundary)
  (when quadtree-boundary
    (second quadtree-boundary)))

(defun quadtree-x1 (quadtree-boundary)
  (when quadtree-boundary
    (third quadtree-boundary)))

(defun quadtree-y1 (quadtree-boundary)
  (when quadtree-boundary
    (fourth quadtree-boundary)))                                               
                             
(defun boundary (quadtree)
  (quadtree-boundary quadtree))

(defun node-p (quadtree)
  (and (quadtree-nw quadtree)
       t))

(defun leaf-p (quadtree)
  (not (node-p quadtree)))

(defun full-p (quadtree)
  (quadtree-object quadtree))

(defun max-depth-p (quadtree)
  (= (quadtree-depth quadtree)
     (quadtree-max-depth quadtree)))

(defun root-p (quadtree)
  (= (quadtree-depth quadtree) 0))

(defun subtrees (quadtree)
  (when quadtree
    (list (quadtree-nw quadtree) (quadtree-ne quadtree)
          (quadtree-sw quadtree) (quadtree-se quadtree))))

(defun intersects-with-subtrees-p (object subtrees)
  (declare (type 3d-vector object))
  (when (and object subtrees)
    ;;(every #'identity (mapcar (lambda (st) (print (point-intersect-p st object)) (point-intersect-p st object)) subtrees))
    (every #'identity (mapcar (lambda (st) (point-intersect-p st object)) subtrees))))

(defun insert-random-object-for-empty (quadtree)
  (when quadtree
    (destructuring-bind (x0 y0 x1 y1) (quadtree-boundary quadtree)
      (setf (quadtree-object quadtree) (make-3d-vector
                                        (/ (abs (- x1 x0)) (1+ (random 10)))
                                        (/ (abs (- y1 y0)) (1+ (random 10)))
                                        (if (full-p quadtree)
                                            (+ *middlepoint-offset* (/ *middlepoint-offset* (1+ (random 10))))
                                            1.0d0))))
    quadtree))

(defun middlepoint-of-quadtree-p (object quadtree)
  (declare (type 3d-vector object))
  (when quadtree
    (if (leaf-p quadtree)
        (destructuring-bind (x0 y0 x1 y1) (quadtree:boundary quadtree)  ;;  d---c
          (let* ((A (make-3d-vector x0 y0 0.0d0))                       ;;  | / |
                 (C (make-3d-vector x1 y1 0.0d0))                       ;;  a---b
                 (M (v+ A
                        (v* ;; calc  middlepoint of boundary
                         (v- c a)
                         0.5d0))))
            (point-intersect-p (make-quadtree (- (x M) *middlepoint-offset*)
                                              (- (y M) *middlepoint-offset*)
                                              (+ (x M) *middlepoint-offset*)
                                              (+ (y M) *middlepoint-offset*))
                               object)))
        (intersects-with-subtrees-p object (subtrees quadtree)))))
     

(defun quadtree-equal (q1 q2)
  ;;TODO (declare (type quadtree q1 q2)))
  
  )
    
(defun subdevide (quadtree &key erease-object)
  (unless (leaf-p quadtree)
    (error "The quadtree ~A is already subdevided." quadtree))
  (destructuring-bind (x0 y0 x1 y1) (boundary quadtree)
    (declare (double-float x0 y0 x1 y1))
    (let ((xmid (/ (+ x1 x0) 2.0))
          (ymid (/ (+ y1 y0) 2.0))
          (depth (1+ (the fixnum (quadtree-depth quadtree))))
          (max-depth (quadtree-max-depth quadtree)))
      (setf (quadtree-nw quadtree) (%make-quadtree-internal x0 ymid xmid y1 depth max-depth)
            (quadtree-ne quadtree) (%make-quadtree-internal xmid ymid x1 y1 depth max-depth)
            (quadtree-sw quadtree) (%make-quadtree-internal x0 y0 xmid ymid depth max-depth)
            (quadtree-se quadtree) (%make-quadtree-internal xmid y0 x1 ymid depth max-depth))))
  (let ((object (quadtree-object quadtree)))
    (if object
        (if erease-object
            (progn
              (insert (quadtree-nw quadtree) (quadtree-empty-middlepoint (quadtree-nw quadtree)))
              (insert (quadtree-ne quadtree) (quadtree-empty-middlepoint (quadtree-ne quadtree)))
              (insert (quadtree-sw quadtree) (quadtree-empty-middlepoint (quadtree-sw quadtree)))
              (insert (quadtree-se quadtree) (quadtree-empty-middlepoint (quadtree-se quadtree)))
              (setf (quadtree-object quadtree) (quadtree-empty-middlepoint quadtree)))
            (progn
              (insert (quadtree-nw quadtree) object)
              (insert (quadtree-ne quadtree) object)
              (insert (quadtree-sw quadtree) object)
              (insert (quadtree-se quadtree) object)))
        (progn
          (insert (quadtree-nw quadtree) (quadtree-empty-middlepoint (quadtree-nw quadtree)))
          (insert (quadtree-ne quadtree) (quadtree-empty-middlepoint (quadtree-ne quadtree)))
          (insert (quadtree-sw quadtree) (quadtree-empty-middlepoint (quadtree-sw quadtree)))
          (insert (quadtree-se quadtree) (quadtree-empty-middlepoint (quadtree-se quadtree)))
          (setf (quadtree-object quadtree) (quadtree-empty-middlepoint quadtree)))))
  quadtree)


(defun coord-intersect-p (quadtree x y)
  (declare (type double-float x y))
  (point-intersect-p quadtree (make-3d-vector x y 0.0d0)))

(defun 2d-dot-product (v-1 v-2)
  "Returns the dot-product"
  (+ (* (x v-1) (x v-2))
     (* (y v-1) (y v-2))))

(defun point-intersect-p (quadtree P)
  (declare (type 3d-vector P))
  (let* ((corner-points (quadtree-boundary->corner-points quadtree))
         (A (first corner-points))  ;;  D---C
         (B (second corner-points)) ;;  | / |
         (C (third corner-points))  ;;  A---B
         (AB (v- B A))
         (AP (v- P A))
         (BC (v- C B))
         (BP (v- P B)))
    (and (<= 0 (2d-dot-product AB AP) (2d-dot-product AB AB))
         (<= 0 (2d-dot-product BC BP) (2d-dot-product BC BC)))))

(defun object-equal (p1 p2)
  (declare (type 3d-vector p1 p2))
  (and (quadtree:same-coords-as p1 p2)
       (same-value-as p1 p2)))

(defun same-coords-as (p1 p2)
  (declare (type 3d-vector p1 p2))
  (and (equal (x p1)
              (x p2))
       (equal (y p1)
              (y p2))))

(defun same-value-as (p1 p2)
  (declare (type 3d-vector p1 p2))
  (equal (z p1)
         (z p2)))
  

(defun set-in-quadtree (quadtree x y)
  (declare (type double-float x y))
  (insert quadtree (make-3d-vector x y 1.0d0)))

(defun get-in-quadtree (quadtree x y &optional neighbor-p)
  (declare (type double-float x y))
  (funcall (lambda (queried)
             (if neighbor-p
                 (when (same-coords-as (make-3d-vector x y 0.0d0) queried)
                   queried)
                 queried))
           (query quadtree x y)))

(defun max-quadtrees (quadtrees)
  "Returns the biggest quadtress which are just the with the biggest
depth value"
  (let ((maximum (alexandria::reduce #'max quadtrees :key #'quadtree-depth)))
    (remove-if-not (alexandria:curry #'equal maximum) quadtrees :key #'quadtree-depth)))

(defun some-subtree-node-p (quadtree)
  "Returns if some subtree from the given `quadtree'
is a node."
  (some #'node-p (subtrees quadtree)))

(defun quadtree-map (quadtree fun &key compose-old-and-new-fun)
  "Maps a function `fun' over the objects saved in `quadtree': meaning
`fun' needs to take one argument which is of type `3d-vector' and
  returns a `3d-vector' too. The input and output of `fun' can then be
  composed with the function `compose-old-and-new-fun' which should a
  `3d-vector'."
  (when (and quadtree fun)
    (loop for subtree in (subtrees quadtree) do
      (with-slots (object) quadtree
        (when (and subtree object)
          (if compose-old-and-new-fun
              (let ((composed (funcall compose-old-and-new-fun
                                       object
                                       (funcall fun object))))
                (declare (type 3d-vector composed))
                (setf object composed))
              (setf object (funcall fun object))))))
    quadtree))

(defun quadtree->list (quadtree &key with-leaves)
  (let ((parents `(,quadtree))
        (visited-parents '()))
    ;; Traverse tree and save visited parents in given var
    (loop while (first parents) do
      (push (pop parents) visited-parents)
      (loop for subtree in (subtrees (first visited-parents)) do
        (when (node-p subtree)
          (setf parents (append parents (list subtree))))))
    (if with-leaves
        (remove-if-not (lambda (st) (and (identity st) (leaf-p st) (full-p st)))
                       (alexandria::reduce #'append (mapcar #'subtrees visited-parents)))
        visited-parents)))

(defun quadtree-leaves (quadtree)
  (remove-if-not #'leaf-p (quadtree->list quadtree)))

(defun minimize-quadtree (quadtree)
  ;; TODO: make quadtree more dynamic by removing blank spaces around
  ;; the middle point by making the quadtree smaller: e. g. from
  ;; 0 0 0 0
  ;; 0 1 1 0
  ;; 0 1 1 0 
  ;; 0 0 0 0
  ;; to
  ;; 1 1
  ;; 1 1
  (let ((nodes (quadtree->list quadtree)))
    ;; Try to remove redundant subtrees
    (loop while nodes do
      (let ((current-max-quadtrees (max-quadtrees nodes)))
        (when (every #'some-subtree-node-p current-max-quadtrees)
          (return))
        (loop for parent in (remove-if #'some-subtree-node-p current-max-quadtrees) do
          (cleanup-leaves parent))
        (setf nodes (nthcdr (length current-max-quadtrees) nodes))))))
;; since quadtree was searched in level order
;; visited-parents is sorted from a high to
;; low depth value

(defun normalize-quadtree (quadtree)
  (let* ((leaves-and-nodes (quadtree->list quadtree :with-leaves t))
         (leaves (remove-if #'node-p leaves-and-nodes))
         (sum 0.0d0))
    (loop for leave in leaves do
      (incf sum (z (quadtree-object leave))))
    (loop for entry in leaves-and-nodes do
      (with-slots (object) entry
        (setf object (make-3d-vector (x object)
                                     (y object)
                                     (float (/ (z object)
                                               sum))))))))


(defun pair-elements (l)
  (when (> (length l) 1)
    (append (list (subseq l 0 2))  (pair-elements (nthcdr 1 l)))))

(defun cleanup-leaves (quadtree)
"If the filled subtrees of `quadtree' have the same value or are all empty,
these subtrees will be removed, since they are redundant."       
  (when quadtree
    (let* ((non-empty-subtrees (remove-if-not #'identity
                                              (subtrees quadtree)
                                              :key #'quadtree-object))
           (object-values-equal (every #'identity (mapcar
                                                   (lambda (l) (same-value-as (first l) (second l)))
                                                   (pair-elements
                                                    (mapcar
                                                    #'quadtree-object
                                                    non-empty-subtrees))))))
      (when (and 
             object-values-equal
             (every #'leaf-p non-empty-subtrees))
        (let ((value (if non-empty-subtrees
                         (quadtree-object (quadtree-nw quadtree))
                         0.0d0)))
          (setf (quadtree-object quadtree)
                (make-3d-vector
                 (x (quadtree-object quadtree))
                 (y (quadtree-object quadtree))
                 (z value)))
          (clear-subtrees-of quadtree))))
    quadtree))


         
(defun clear-in-quadtree (quadtree x y &optional remove-subtrees-if-possible)
  ;; Removes point with same coords and if subtrees of point are empty,
  ;; the the subtrees will be removed to save space: therfore the parent will
  ;; be updated
  (declare (type double-float x y))
  (let ((cleared-point (make-3d-vector x y 0.0d0)))
    (multiple-value-bind (queried-point queried-trees) (query quadtree x y t)
      (when (same-coords-as cleared-point queried-point)
        (setf (quadtree-object (first (last queried-trees 1))) cleared-point)
        (when remove-subtrees-if-possible
          (let ((parent (first (last queried-trees 2))))
            (remove-redundant-subtrees parent)))))))

(defun round-to (number precision &optional (what #'round))
    (let ((div (expt 10 precision)))
      (/ (funcall what (* number div)) div)))

(defun add-offset-if-middlepoint (object quadtree)
  (declare (type 3d-vector object))
  (when (middlepoint-of-quadtree-p object quadtree)
    (setf object (make-3d-vector (round-to (- (x object) *middlepoint-offset*) 5)
                                 (round-to (- (y object) *middlepoint-offset*) 5)
                                 (z object))))
  object)

(defun quadtree-empty-middlepoint (quadtree)
  (when quadtree
    (destructuring-bind (x0 y0 x1 y1) (boundary quadtree)      
      (let ((START (make-3d-vector x0 y0 0.0d0))
            (END (make-3d-vector x1 y1 0.0d0)))
        (v+ START
            (v*
             (v- END START)
             0.5d0))))))


(defun insert (quadtree object &key wrt-resolution)
  (cond
    ;; When the object does not intersect the quadtree, just return nil.
    ((not (if quadtree (point-intersect-p quadtree object))) nil)
    ;; When the quadtree is a node, recursively insert the object to its children.
    ((node-p quadtree)
     (insert (quadtree-nw quadtree) object :wrt-resolution wrt-resolution)
     (insert (quadtree-ne quadtree) object :wrt-resolution wrt-resolution)
     (insert (quadtree-sw quadtree) object :wrt-resolution wrt-resolution)
     (insert (quadtree-se quadtree) object :wrt-resolution wrt-resolution)
     t)
    ;; Insert the object to the quadtree, if leaf is free or if the
    ;; resolution of the quadtree is smaller or equal to the given
    ;; resolution.
    ((and (not (full-p quadtree))
          (if wrt-resolution
              (<= (quadtree-resolution quadtree) wrt-resolution)
              t))
     (setf (quadtree-object quadtree) (add-offset-if-middlepoint object quadtree))
     t)
    ;; If leaf is not free, but has same coords as object or
    ;; resolution of this quadtree is smaller or equal to the given
    ;; resolution, update the object
    ((and (full-p quadtree)
          (or
           (if wrt-resolution
               (<= (quadtree-resolution quadtree) wrt-resolution)
               nil)
           (same-coords-as (quadtree-object quadtree) object)))
     (setf (quadtree-object quadtree) object)
     t)
    ;; When the quadtree is full and is not at its max depth, subdevide it and
    ;; recursively insert the object.
    ((and (not (max-depth-p quadtree))
               (if wrt-resolution
                  (> (quadtree-resolution quadtree) wrt-resolution)
                  t))         
     (subdevide quadtree :erease-object wrt-resolution)
     (insert quadtree object :wrt-resolution wrt-resolution))
    ;; Otherwise if the max-depth-p or the objects were equal (object-equal returned t),
    ;; do nothing, since the object represents accordingly the given object. 
    (t
     t)))

(defun v-abs (v)
  (declare (type 3d-vector v))
  (make-3d-vector (abs (x v)) (abs (y v)) (abs (z v))))

(defun rotate-quadtree (quadtree angle &optional M-root)
  (destructuring-bind (x0 y0 x1 y1) (boundary quadtree)
    (flet ((rotate-v (v)
             (rotate 
              (axis-angle->quaternion 
               (make-3d-vector 0 0 1)
               angle)
              v)))
      (let* ((START (make-3d-vector x0 y0 0.0d0))
             (END (make-3d-vector x1 y1 0.0d0))
             (M (v-abs (if M-root
                           M-root
                           (v+
                            START
                            (v*
                             (v- END START)
                             0.5d0)))))
             (M-START (v- START M))
             (M-END (v- END M))
             (rotated-start (v+ (rotate-v M-START)
                                M))
             (rotated-end (v+ (rotate-v M-END)
                              M)))
        (print M-START)
        (print M-END)
        (setf (quadtree-boundary quadtree)
              (list (x rotated-start)
                    (y rotated-start)
                    (x rotated-end)
                    (y rotated-end)))
        (when (quadtree-object quadtree) ;; TODO: check if below works 100%
          (with-slots (object) quadtree
            (let* ((M-of-quadtree (v*
                                   (v- END START)
                                   0.5d0))
                   (M-TO-M-of-quadtree (v- M M-of-quadtree))
                   (rotated-obj (v+
                                 START
                                 M
                                 (rotate-v M-TO-M-of-quadtree))))
              (setf object (make-3d-vector (x rotated-obj) (y rotated-obj) (z object))))))
        (mapcar (alexandria:rcurry #'rotate-quadtree angle M)
                (remove-if-not #'identity (subtrees quadtree)))))))
        
    
  


(defun query (quadtree x y &optional trees)
  (declare (type double-float x y))
  (let ((ret (query-tree quadtree x y '())))
    (when ret
      (values (quadtree-object (first (last ret)))
              (when trees
                ret)))))

(defun query-tree (quadtree x y path)
  (declare (type double-float x y))
  (when (coord-intersect-p quadtree x y)
    (if (node-p quadtree)
        (or (query-tree (quadtree-nw quadtree) x y (append path (list (quadtree-nw quadtree))))
            (query-tree (quadtree-ne quadtree) x y (append path (list (quadtree-ne quadtree))))
            (query-tree (quadtree-sw quadtree) x y (append path (list (quadtree-sw quadtree))))
            (query-tree (quadtree-se quadtree) x y (append path (list (quadtree-se quadtree)))))
        (when (quadtree-object quadtree)
          path))))

(defun clear-subtrees-of (quadtree)
  (setf (quadtree-nw quadtree) nil
        (quadtree-ne quadtree) nil
        (quadtree-sw quadtree) nil
        (quadtree-se quadtree) nil)
  t)

(defun clear-quadtree (quadtree)
  (setf (quadtree-object quadtree) nil)
  (clear-subtrees-of quadtree)
  t)

(defun matrix->quadtree (origin-x origin-y resolution matrix &optional (threshold 0.0d0))
  "Creates an quadtree from the matrix `matrix'. The elements from `matrix' will
be mapped to the centers of the points. Since adding zeros is important too, this
will be done too. `origin-x' and `origin-y' have to be the bottom left corner and
`resolution' is the width of a cell in the returned quadtree."
  (declare (type double-float origin-x origin-y resolution)
           (type (simple-array * 2) matrix))
  (let* ((width (coerce (* (array-dimension matrix 1) resolution) 'double-float))
         (height (coerce (* (array-dimension matrix 0) resolution) 'double-float))
         (x0 origin-x)
         (x1 (+ origin-x width))
         (y0 origin-y)
         (y1 (+ origin-y height))
         (qt (make-quadtree x0 y0 x1 y1)))
    (dotimes (y (array-dimension matrix 1))
      (let ((backward-y (- (1- (array-dimension matrix 1)) y)))
        (dotimes (x (array-dimension matrix 0))
          (when (>= (aref matrix backward-y x) threshold)
            (insert qt
                    (make-3d-vector
                     (+ (/ resolution 2) 
                        (+ origin-x (* x resolution)))
                     (- (+ origin-y (* y resolution) resolution)
                        (/ resolution 2))
                     (aref matrix backward-y x)))))))
    qt))

(defun base-nth (l1 l2 n fun)
  (if (funcall fun (nth n l1) (nth n l2))
      (nth n l1)
      (nth n l2)))

(defun >nth (l1 l2 n)
  (base-nth l1 l2 n #'>))

(defun <nth (l1 l2 n)
  (base-nth l1 l2 n #'<))

(defun quadtree-resolution (qt)
  (destructuring-bind (x0 y0 x1 y1) (quadtree-boundary qt)
    (let ((res0 (round-to (abs (- x1 x0)) 5))
          (res1 (round-to (abs (- y1 y0)) 5)))
      (if (equal res0 res1)
          res0
          (error 'simple-error "oh no")))))

(defun merge-quadtrees (qt1 &rest quadtrees)
  "TODO"
  (if (first quadtrees)
      (first quadtrees)
      qt1))
  ;; (when quadtrees                       
  ;;   (let* ((merged-qt '())
  ;;          (leaves-qt1 (quadtree-leaves qt1))
  ;;          (leaves-qt2 (quadtree-leaves (first quadtrees)))
   
  ;;          ))))

(defun quadtree-boundary->vectors (qt)
  "Returns the boundary of the given quadtree `qt' as vectors AB, BC,
CD, DA in a list."
  (mapcar #'v- (mapcar #'reverse (quadtree-boundary->edge-points qt))))

(defun quadtree-boundary->edge-points (qt)
  "Returns the edges of the quadtree `qt' in a list containing lists
  with start and end points (A, B, C or D) from each edge."
  (when qt
    (let* ((corner-points (quadtree-boundary->corner-points qt))
           (A (first corner-points))
           (B (second corner-points))
           (C (third corner-points))
           (D (fourth corner-points)))
      (list `(,A ,B) `(,B ,C) `(,C ,D) `(,D ,A)))))

(defun quadtree-boundary->corner-points (qt)
  "Returns the points of the corners `A', `B', `C' and `D' of given
quadtree `qt' in a list."
  (when qt
    (destructuring-bind (x0 y0 x1 y1) (quadtree-boundary qt)        ;;  D---C
      (let* ((A (make-3d-vector x0 y0 0.0d0))                       ;;  | M |
             (C (make-3d-vector x1 y1 0.0d0))                       ;;  A---B
             (AM (v* ;; calc  middlepoint of boundary
                  (v- C A)
                  0.5d0))
             (B (v+
                 A AM
                 (rotate ;; calc point M, which is between a and c 
                  (axis-angle->quaternion 
                   (make-3d-vector 0 0 1)
                   (+ pi (/ pi 2))) ;; rotate around 270°C
                  AM)))
             (D (v+
                 A AM
                 (rotate ;; calc point M, which is between a and c 
                  (axis-angle->quaternion 
                   (make-3d-vector 0 0 1)
                   (/ pi 2)) ;; rotate around 90°C
                  AM))))
        (mapcar #'v-round-to (list A B C D))))))

(defun points->vector-points (points)
  "Returns all (redundant) vectors/connections between each point in
  `points'. The vectors are encoded as a list of a start and
   end point. Therefore, a list of lists will be returned."
  (when points
    (let ((ret '()))
      (loop for i from 0 to (1- (length points)) do
        (loop for j from 0 to (1- (length points)) do
          (when (not (equal j i))
            (push (list (nth i points) (nth j points)) ret))))
      ret)))

(defun intersecting-vector-points-lists (vec-points-list-1
  vec-points-list-2 &optional intersecting-counts-only-with-limited-vector-length-p)
  "Gets two lists `vec-points-list-1' and `vec-points-list-2' which
  are lists containing vectors in form of a start and end point. Meaning
  vector AB is here a list: (list (point a) (point b)). Between each
  of these sublists in `vec-points-list-1' and `vec-points-list-2' a
  intersection point will be calculated (or the given
  sublists/vectors are parallel to each other). A list of these
  intersection points will be returned. If
  `intersecting-counts-only-with-limited-vector-length-p' is T the
  intersection between the given vectors/sublists has to be valid."
  (when (and vec-points-list-1 vec-points-list-2)
    (let ((ret '()))
      (loop for vec-point-1 in vec-points-list-1 do
        (loop for vec-point-2 in vec-points-list-2 do
          (multiple-value-bind (intersected-point valid-intersection-p)
              (point-of-intersection (first vec-point-1) (second vec-point-1)
                                     (first vec-point-2) (second vec-point-2))
            (when intersected-point
              (let ((round-intersected-point (v-round-to intersected-point)))
                (when (or (not ret)
                          (not (member round-intersected-point ret :test #'v-equal)))
                  (if intersecting-counts-only-with-limited-vector-length-p
                      (when valid-intersection-p
                        (push round-intersected-point ret))
                      (push round-intersected-point ret))))))))
      ret)))

(defun max-points-on-vector-points (ps vec-points)
  "Takes points from `ps' and returns the points `min' and `max' in a
  list, which are the farthest away from each other on the given vector
  `vec-points'. `vec-points' is represented as a list of start (first)
  and end (second) point of the vector."
  (when (and ps (equal 2 (length vec-points)))
    (let* ((points-on-vec (mapcar
                          #'v-round-to
                          (remove-if-not
                           (alexandria:rcurry #'point-on-vector-points-p vec-points)
                           ps)))
          (min (first points-on-vec))
          (max (first points-on-vec)))
      (when points-on-vec
      (loop for p in points-on-vec do
        (when (and (<= (x p) (x min)) 
                   (<= (y p) (y min)))
          (setf min p))
        (when (and (>= (x p) (x max))
                   (>= (y p) (y max)))
          (setf max p)))
      (list min max)))))

(defun point-on-vector-points-p (p vec-points)
  "Returns T if given point `p' is on the vector `vec-points' which is
  represented as a list of a start (first) and end point (second) of the vector."
  (declare (type 3d-vector p))
  (when (equal 2 (length vec-points))
    (point-on-vector-p (v- p (first vec-points))
                       (v- (second vec-points) (first vec-points)))))
  
(defun point-on-vector-p (p vec)
  "Returns T if given point `p' is on the given vector `vec'."
  (declare (type 3d-vector p vec))
  (cond
    ;; if (x vec) is zero
    ((and (equal (x vec) 0.0d0)
          (not (equal (y vec) 0.0d0)))
     (and (equal (x p) 0.0d0)
          (<= 0.0d0 (/ (y p) (y vec)) 1.0d0)))
    ;; if (y vec) is zero
    ((and (equal (y vec) 0.0d0)
          (not (equal (x vec) 0.0d0)))
     (and (equal (y p) 0.0d0)
          (<= 0.0d0 (/ (x p) (x vec)) 1.0d0)))
    ;; if (x vec) and (y vec) are zero
    ((and (equal (x vec) 0.0d0)
          (equal (y vec) 0.0d0))
     (and (equal (x p) 0.0d0)
          (equal (y p) 0.0d0)))
    ;; if neither (x vec) nor (y vec) are zero
    (t
     (let ((s (/ (x p) (x vec)))
           (r (/ (y p) (y vec))))
       (and (equal s r)
            (<= 0.0d0 s 1.0d0))))))

(defun intersecting-points (qt1 qt2)
  "Calculates the points bordering the area which the quadtrees
  `qt1' and `qt2' both cover."
  (when (and qt1 qt2
             (leaf-p qt1)
             (leaf-p qt2))
    (let* ((qt1-edge-points (quadtree-boundary->edge-points qt1))
           (qt2-edge-points (quadtree-boundary->edge-points qt2))
           ;; calculates the intersection points between the edges
           ;; of qt1 and qt2
           (intersecting-points-between-edges
             (intersecting-vector-points-lists qt1-edge-points
                                               qt2-edge-points))
           ;; creates new vectors from the intersection points from
           ;; the edges of qt1 and qt2 <- bull
           ;;
           ;; TODO: BETTER: ADD corners if in both qts <- only smaller
           ;; bull 
           (vector-points-of-intersecting-points
             (points->vector-points
              (append intersecting-points-between-edges
                      (quadtree-boundary->corner-points qt1)
                      (quadtree-boundary->corner-points qt2))))
           ;; calculates all points which intersect with the edges of
           ;; qt1 and qt2 
           (all-intersecting-points
             (append
              (intersecting-vector-points-lists
               vector-points-of-intersecting-points
               qt2-edge-points
               t)
              (intersecting-vector-points-lists
               vector-points-of-intersecting-points
               qt1-edge-points
               t)))
           ;; removes all points which do not intersect within the
           ;; boundaries of qt1 and qt2
           (intersecting-points-on-edges
             (remove-if-not (lambda (p)
                              (and (point-intersect-p qt1 p)
                                   (point-intersect-p qt2 p)))                           
                            all-intersecting-points))
           ;; calculates the points which are the farthest away
           ;; on each edge of qt1 or qt2
           (max-intersecting-points-on-edges
             (remove-if-not
              #'identity
              (mapcar (alexandria:curry
                       #'max-points-on-vector-points
                       intersecting-points-on-edges)
                      (append qt1-edge-points
                              qt2-edge-points)))))
      (let ((ret '()))
        (loop for pair in max-intersecting-points-on-edges
              collect (loop for p in pair
                            collect (unless (member p ret :test #'v-equal)
                                      (push p ret))))
        ret))))


(defun v-equal (v w)
  "Returns T, if x, y and z are equal between vector `v' and `w'."
  (declare (type 3d-vector v w))
  (and (equal (x v) (x w))
       (equal (y v) (y w))
       (equal (z v) (z w))))

(defun v-round-to (v &optional (precision 3))
  "Rounds x and y of given vetor `v' w.r.t. the `precision'."
  (declare (type 3d-vector v))
  (make-3d-vector (round-to (x v) precision)
                  (round-to (y v) precision)
                  (z v)))

(defun point-of-intersection (a b c d)
  "This function calculates the possible intersection point between the vector AB
  and the vector CD: meaning `a' is the start point of AB and `b' is the
  end point (analog for `c', `d' and CD). It returns the intersecting
  point and if the vectors intersected or the extension of
  these. Explanation of the calculation is here: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C6F1EC3B8FA75FB4FADC3DBA35DE4D6C/9780511804120c7_p220-293_CBO.pdf/search_and_intersection.pdf"
  (declare (type 3d-vector a b c d))
  (unless (or (v-equal a c)
              (v-equal b d))
    (let ((denom (+ (* (x a) (- (y d) (y c)))
                    (* (x b) (- (y c) (y d)))
                    (* (x d) (- (y b) (y a)))
                    (* (x c) (- (y a) (y b))))))
      (unless (equal denom 0.0d0)
        (let* ((num-s (+ (* (x a) (- (y d) (y c)))
                         (* (x c) (- (y a) (y d)))
                         (* (x d) (- (y c) (y a)))))
               (num-r (* -1 
                         (+ (* (x a) (- (y c) (y b)))
                            (* (x b) (- (y a) (y c)))
                            (* (x c) (- (y b) (y a))))))
               (s (/ num-s denom))
               (r (/ num-r denom)))
         
          (values
           (make-3d-vector (+ (x a) 
                              (* s
                                 (- (x b) (x a))))
                           (+ (y d) 
                              (* s
                                 (- (y b) (y a))))
                           1.0d0)
           (and (<= 0.0d0 s 1.0d0)
                (<= 0.0d0 r 1.0d0))))))))
