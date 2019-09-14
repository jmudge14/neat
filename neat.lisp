;;;; neat.lisp

(in-package #:neat)

;;;; Utilities

(defmacro defincrementer (name)
  (let ((var (intern (string-upcase 
                       (concatenate 'string "*" (symbol-name name) "*"))))
        (func (intern (string-upcase
                        (concatenate 'string "next-" (symbol-name name))))))
    `(progn (defvar ,var 0)
            (defun ,func ()
              (setf ,var (1+ ,var))))))

(defmacro extendf (place object)
  "Adds object to the end of list via append"
  `(setf ,place (append ,place (list ,object))))


;;;; NEAT implementation

; Global innovation number
(defincrementer innovation)

; Generation numbers - for aging genomes
(defincrementer generation)

; Node ID's - for assigning new nodes somewhat sensibly
(defincrementer node-id)


; Node reperesentation
(defclass node ()
  ((id :initarg :id :reader node-id)
   (type :initarg :type :reader node-type
         :documentation "One of :sensor :output :bias or :hidden.")
   (value :initform 0.0 :accessor node-value
          :documentation "Current value of the node. For hidden nodes, this is the intermediate
                          sum during evaluation. For all others, this is the constant value they
                          will output.")))

(defvar *nodes* nil)
(defun make-node (type &optional (node-id (next-node-id)))
  (let ((node (make-instance 'node
                             :id node-id
                             :type type)))
    (push node *nodes*)))


(defmethod print-object ((node node) stream)
  (format stream "#<NODE ~A:~A>" 
          (symbol-name (node-type node))
          (node-id node)))

(defclass connection ()
  ((innovation :initarg :innovation
               :reader innovation)
   (in-node :initarg :in-node :reader in-node)
   (out-node :initarg :out-node :reader out-node)
   (weight :initarg :weight :accessor weight)
   (enabled :initarg :enabled :accessor enabled)))

(defmethod print-object ((connection connection) stream)
  (with-slots (innovation in-node out-node weight enabled) connection
    (format stream "#<CONN ~A ~A~A->~A @~A>"
            innovation
            (if enabled "" "DIS ")
            (node-id in-node)
            (node-id out-node)
            weight)))

(defun make-connection (in-node out-node weight 
                        &key (enabled t) 
                             (innovation (next-innovation)))
  (make-instance 'connection
                 :innovation innovation
                 :in-node in-node
                 :out-node out-node
                 :weight weight
                 :enabled enabled))

(defmethod copy-connection ((connection connection))
  (let ((in-node (in-node connection))
        (out-node (out-node connection))
        (weight (weight connection))
        (enabled (enabled connection))
        (innovation (innovation connection)))
    (make-instance 'connection :in-node in-node
                               :out-node out-node
                               :weight weight
                               :enabled enabled
                               :innovation innovation)))

(defclass genome ()
  ((connections :accessor connections :initarg :connections :initform nil)
   (nodes :accessor nodes :initarg :nodes :initform nil)
   (generation :initarg :generation :reader generation :initform 0)
   (fitness :initform 0.0 :initarg :fitness :accessor fitness)))

(defmethod copy-genome ((genome genome))
  "Copy a genome and set its generation number to current"
  (let ((connections (map 'list #'copy-connection (connections genome)))
        (nodes (nodes genome)))
    (make-instance 'genome :generation *generation* 
                           :connections connections
                           :nodes nodes)))

(defparameter +connection-perturb-probability+ 0.90
  "Chance of any particular connection being perturbed")
(defparameter +connection-perturb-scale+ 0.10
  "Scale factor of weight against normal distribution")

(defmethod perturb ((genome genome))
  "Randomly adjust connection weights"
  (dolist (connection (connections genome))
    (let* ((new-weight (+ (weight connection)
                          (* +connection-perturb-scale+
                             (alexandria:gaussian-random))))
           (perturb? (> +connection-perturb-probability+ 
                        (random 1.0))))
      (when perturb?
        (setf (weight connection) new-weight)))))

(defmacro where (&rest l) `(lambda (it) (or ,@l)))

(defmethod find-connection ((genome genome) (n1 node) (n2 node))
  (find (list n1 n2) 
        (connections genome)
        :test #'equalp
        :key (lambda (connection) 
               (list (in-node connection)
                     (out-node connection)))))


(defmethod find-new-connection ((genome genome))
  "Return a pair of nodes that are not currently connected. 
   If such a pair can't be found after too many tries, return nil."
  (loop for i from 0
        for n1 = (alexandria:random-elt (nodes genome))
        and n2 = (alexandria:random-elt (nodes genome))
        until (or (> i 100)
                  (and (not (equalp n1 n2))
                       (not (find-connection genome n1 n2))
                       ; no connections *to* an input node (from is allowed)
                       (not (find (node-type n2) '(:sensor :bias)))))
        finally (if (<= i 99)
                    (return (list n1 n2))
                    (return nil))))

(defparameter +new-connection-weight-magnitude+ 1.0
  "Largest absolute value of new connections' random weights")

(defun signed-random (magnitude)
  (let ((sign (> (random 1.0) 0.5)))
    (if sign
        (- (random magnitude))
        (random magnitude))))


(defmethod add-connection ((genome genome) &optional nodes weight)
  "Add connection mutation - new connection gene with random weight is added connecting 
   two previously unconnected nodes. The optional nodes and weight allow this to be used
   during minimal genome construction."
  (let ((new-nodes (or nodes 
                       (find-new-connection genome))))
    (when new-nodes
      (let* ((in-node (first new-nodes))
             (out-node (second new-nodes))
             (new-connection 
               (make-connection in-node
                                out-node
                                (or weight 
                                    (signed-random +new-connection-weight-magnitude+)))))
        (extendf (connections genome) new-connection)))))

(defmethod add-node ((genome genome))
  "Add node mutation - existing connection is split and new node placed where old connection 
   used to be. Old connection is disabled and two new connections are added to the genome."
  (let* ((old-conn (alexandria:random-elt (connections genome)))
         ; nodes
         (in-node (in-node old-conn))
         (new-node (make-node :hidden))
         (out-node (out-node old-conn))
         ; new connections
         (in-conn (make-connection in-node new-node 1.0))
         (out-conn (make-connection new-node out-node (weight old-conn))))
    (setf (enabled old-conn) nil)
    (extendf (nodes genome) new-node)
    (extendf (connections genome) in-conn)
    (extendf (connections genome) out-conn)))

(defmethod add-node-of-type ((genome genome) node-type &optional (node-id (next-node-id)))
  (let* ((new-node (make-node node-type node-id))
         (split-pos (position-if (alexandria:rcurry #'eq node-type)
                                 (nodes genome)
                                 :key #'node-type
                                 :from-end t))
         (head (if split-pos
                   (subseq (nodes genome) 0 (1+ split-pos))
                   (nodes genome)))
         (tail (if split-pos
                   (subseq (nodes genome) (1+ split-pos))
                   nil)))
    (setf (nodes genome)
          (append head (list new-node) tail))))


;;;; Breeding routines

(defmethod match-genes ((gen-a genome) (gen-b genome))
  (do* ((conns-a (connections gen-a))
        (conns-b (connections gen-b))
        (matches nil)
        (count-match 0)
        (count-disjoint 0)
        (count-excess 0)
        (conn-a (pop conns-a))
        (conn-b (pop conns-b)))
      ((not (or conn-a conn-b)) ; Stop when we have exhausted all genes
       (values (nreverse matches) 
               count-match 
               count-disjoint 
               count-excess))
      (let ((ia (and conn-a (innovation conn-a)))
            (ib (and conn-b (innovation conn-b))))
        (cond ((and ia ib (= ia ib))
               (push (cons conn-a conn-b) matches)
               (incf count-match)
               (setf conn-a (pop conns-a)
                     conn-b (pop conns-b)))
              ((and ia ib (< ia ib))
               (push (cons conn-a nil) matches)
               (incf count-disjoint)
               (setf conn-a (pop conns-a)))
              ((and ia ib (< ib ia))
               (push (cons nil conn-b) matches)
               (incf count-disjoint)
               (setf conn-b (pop conns-b)))
              ((or (not ia)
                   (not ib))
               (push (cons conn-a conn-b) matches) ; one will be nil
               (incf count-excess)
               (setf conn-a (pop conns-a)
                     conn-b (pop conns-b)))))))


(defparameter +distance-excess-coefficient+ 1.0)
(defparameter +distance-disjoint-coefficient+ 1.0)
(defparameter +distance-weights-coefficient+ 1.0)

(defmethod distance ((gen-a genome) (gen-b genome))
  (multiple-value-bind (matches count-match count-disjoint count-excess)
                       (match-genes gen-a gen-b)
    (labels ((weight-diff (match)
               (let ((m1 (car match))
                     (m2 (cdr match)))
                 (or (and m1 m2 (abs (- (weight m1)
                                        (weight m2))))
                     0))))
      (let* ((matched-weight-differences (reduce #'+ (map 'list #'weight-diff matches)))
             (average-weight-differences (/ matched-weight-differences count-match))
             (num-genes (max (length (connections gen-a))
                             (length (connections gen-b))
                             1)))
        (+ (/ (* count-excess +distance-excess-coefficient+) num-genes)
           (/ (* count-disjoint +distance-disjoint-coefficient+) num-genes)
           (* average-weight-differences +distance-weights-coefficient+))))))

(defparameter +compatibility-distance-threshold+ 1.0)



(defparameter +equal-fitness-threshold+ 0.2)

(defmethod breed ((gen-a genome) (gen-b genome))
  "Return a new genome by breeding gen-a and gen-b."
  (let* ((matches (match-genes gen-a gen-b))
         (selector-func (cond ((< (abs (- (fitness gen-a) ; If 'equal' fitness, take a parent randomly.
                                          (fitness gen-b)))
                                  +equal-fitness-threshold+)
                               (if (< (random 1.0) 0.5)
                                   #'car #'cdr))
                              ((< (fitness gen-a)
                                  (fitness gen-b))
                               #'car)
                              (t #'cdr)))
         ; Genome construction
         (connections (map 'list #'copy-connection 
                           (remove-if #'null
                                 (map 'list selector-func
                                      matches))))
         (nodes nil))
    ; Identify all nodes
    (dolist (conn connections)
      (pushnew (in-node conn) nodes)
      (pushnew (out-node conn) nodes))
    ; Yay!
    (make-instance 'genome
                   :connections connections
                   :nodes nodes
                   :generation *generation*)))

(defclass species ()
  ((population :initform nil :initarg :population :accessor population)
   (representative :initform nil :initarg :representative :accessor representative)))

(defmethod make-species ((genome genome))
  (make-instance 'species
                 :population (list genome)
                 :representative genome))

(defmethod speciesp ((species species) (genome genome))
  (< (distance genome
               (representative species))))


(defclass generation ()
  ((species :initform nil :initarg :species :accessor species
            :documentation "List of species present in this population.")
   (generation-id :initform (next-generation) :accessor generation-id)))



