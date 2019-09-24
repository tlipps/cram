#|
  This file is a part of quadtree project.
  Copyright (c) 2015 Masayuki Takagi (kamonama@gmail.com)
|#

(in-package :cl-user)
(defpackage quadtree
  (:use :cl :cl-transforms)
  (:export #:make-quadtree-internal
           #:make-quadtree
           #:boundary
           #:insert
           #:query
           #:clear-quadtree
           #:point-intersect-p
           #:object-equal
           #:same-coords-as))
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
    (print subtrees)
    (print object)
    (every #'identity (mapcar (lambda (st) (print (point-intersect-p st object)) (point-intersect-p st object)) subtrees))
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
    
(defun subdevide (quadtree)
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
        (progn
          (insert (quadtree-nw quadtree) object)
          (insert (quadtree-ne quadtree) object)
          (insert (quadtree-sw quadtree) object)
          (insert (quadtree-se quadtree) object))
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
  (destructuring-bind (x0 y0 x1 y1) (quadtree:boundary quadtree)  ;;  d---c
    (let* ((A (make-3d-vector x0 y0 0.0d0))                       ;;  | / |
           (C (make-3d-vector x1 y1 0.0d0))                       ;;  a---b
           (AM (v* ;; calc  middlepoint of boundary
                (v- c a)
                0.5d0))
           (B (v+
               A AM
               (rotate ;; calc point b, which is between a and c 
                (axis-angle->quaternion 
                 (make-3d-vector 0 0 1)
                 (+ pi (/ pi 2))) ;; rotate around 270°C
                AM)))
           (AB (v- B A))
           (AP (v- P A))
           (BC (v- C B))
           (BP (v- P B)))
      (and (<= 0 (2d-dot-product AB AP) (2d-dot-product AB AB))
           (<= 0 (2d-dot-product BC BP) (2d-dot-product BC BC))))))

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
    ;; If leaf is not free, but has same coords as object, update the object
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
     (subdevide quadtree)
     (print "after subdiv")
     (print quadtree)
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
             (M (v-abs (if M-root ;; M is actually the vector from START to M-root
                           M-root
                           (v*
                            (v- END START)
                            0.5d0))))
             (M-START (v- START M))
             (M-END (v- END M))
             (rotated-start (v+ (rotate-v M-START)
                                M))
             (rotated-end (v+ (rotate-v M-END)
                              M)))
        (setf (quadtree-boundary quadtree)
              (list (x rotated-start)
                    (y rotated-start)
                    (x rotated-end)
                    (y rotated-end)))
        (when (quadtree-object quadtree)
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
              (if trees
                  ret
                  nil)))))

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
          (error #'simple-error "oh no")))))

(defun merge-quadtrees (qt1 &rest quadtrees)
  (let* ((sorted-by-resolution (sort (append `(,qt1) quadtrees) #'<=
                                     :key #'quadtree-resolution))
         (boundries (mapcar #'quadtree-boundary (append `(,qt1)
                                                        quadtrees)))
         (min (reduce (lambda (l1 l2)
                        (list (<nth l1 l2 0)
                              (<nth l1 l2 1)))
                      boundries))
         (max (reduce (lambda (l1 l2)
                        (list (>nth l1 l2 2)
                              (>nth l1 l2 3)))
                      boundries))
         (merged-qt (make-quadtree (first min) (second min)
                                   (first max) (second max))))
   
  ))
