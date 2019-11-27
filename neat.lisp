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
         :documentation "One of :sensor :output :bias or :hidden.")))

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
   (fitness :initform 0.0 :initarg :fitness :accessor fitness)
   (network :initform nil :accessor network)
   (last-update :initform -1 :accessor last-update)))

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
  (dolist (connection (connections genome) genome)
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
        (extendf (connections genome) new-connection)))
    genome))

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
    (extendf (connections genome) out-conn)
    genome))

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
      (let* ((matched-weight-differences (sum (map 'list #'weight-diff matches)))
             (average-weight-differences (/ matched-weight-differences count-match))
             (num-genes (max (length (connections gen-a))
                             (length (connections gen-b))
                             1)))
        (+ (/ (* count-excess +distance-excess-coefficient+) num-genes)
           (/ (* count-disjoint +distance-disjoint-coefficient+) num-genes)
           (* average-weight-differences +distance-weights-coefficient+))))))


(defun sum (list)
  (reduce #'+ list :initial-value 0))


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

(defparameter +compatibility-distance-threshold+ 1.0 
  "Maximum distance for two genomes to be of the same species.")

(defmethod speciesp ((species species) (genome genome))
  (< (distance genome
               (representative species))
     +compatibility-distance-threshold+))


(defmethod initialize ((species species))
  "Return a species with a random representative"
  (make-instance 'species 
                 :population nil
                 :representative (alexandria:random-elt (population species))))

(defmethod copy-species ((species species))
  (make-instance 'species
                 :population (mapcar #'copy-genome (population species))
                 ; Note - representative is read-only in all other cases
                 ; so it's not being copied, it only matters for distance calculations.
                 :representative (representative species)))

(defun speciate (species-list population)
  "Place each member of the population into their respective species, or a new species.
   Return a list of (possibly modified) species."
  (loop for genome in population do
        (unless (loop for species in species-list
                      when (speciesp species genome) 
                        do (push genome (population species))
                        and return t)
          (setf species-list
                (nconc species-list
                       (list (make-instance 'species
                                            :population (list genome)
                                            :representative genome))))))
  species-list)

(defclass generation ()
  ((species :initform nil :initarg :species :accessor species
            :documentation "List of species present in this population.")
   (generation-id :initform (next-generation) :accessor generation-id)
   (population-size :initform 150 :initarg :population-size :accessor population-size)))


(defun make-generation (genome-list &optional (species-list nil) (population-size 150))
  "Make a generation out of the explicit population genome-list.
   When speciating, use representatives from species-list."
  (let* ((population (mapcar #'copy-genome genome-list))
         (species (mapcar #'initialize species-list)))
    (setf species (speciate species population))
    (make-instance 'generation 
                   :species species
                   :population-size population-size)))

(defun weighted-random-elt (weighted-alist)
  "Select a key value at random from weighted-alist.
   Weights are relative - the sum of them is always 100%."
  (let* ((total (sum (mapcar #'cdr weighted-alist)))
         (selected (random total)))
    (dolist (item weighted-alist nil)
      (if (<= selected (cdr item))
          (return-from weighted-random-elt (car item))
          (setf selected (- selected (cdr item)))))))

(defun do-until (predicate func &rest args)
  "Continually apply args to func until the predicate is true against the result,
   and return the result for which predicate was true."
  (loop for r = (apply func args) then (apply func args)
        when (funcall predicate r) return r))

(defun select-random-n-elt (list n)
  "Select n distinct elements from the given list."
  (assert (>= (length list) n) (list n) "Too few elements in list")
  (let ((result nil)
        (remaining list))
    (labels ((not-in-result-list (element)
               (not (find element result))))
      (dotimes (i n result)
        (let ((next-choice (do-until #'not-in-result-list #'alexandria:random-elt remaining)))
          (setf remaining (remove next-choice remaining))
          (push next-choice result))))))

(defparameter +species-mutation-rates+
  '((breed          . 100)  ; (breed g1 g2)
    (weight         . 100)  ; (perturb g1)
    (new-node       . 100)  ; (add-node g1)
    (new-connection . 100)  ; (add-connection g1)
    (copy           . 100)) ; (copy-genome g1)
  "Relative weight of various mutation types. 
   Copy does not mutate. Breed is sexual reproduction.
   All others are asexual reproduction.")

(defmethod reproduce ((species species) &rest args)
  "Return count new genomes from the given species"
  (let* ((count (or (first args) 0)) ; Function parameter - number of new organisms to generate 
         (new-population nil))
    ; Create count new genomes
    (dotimes (c count)
      (let* ((action (if (<= 1 (length (population species)))
                         ; When we only have one specimen, we can't breed, select one of the others.
                         (do-until (lambda (a) (not (eq a 'breed))) 
                                   #'weighted-random-elt 
                                   +species-mutation-rates+)
                         (weighted-random-elt +species-mutation-rates+)))
             (genomes (select-random-n-elt (population species) 
                                           (if (eq action 'breed) 2 1))))
        (push (case action
                (breed          (breed (first genomes) (second genomes)))
                (weight         (perturb (copy-genome (first genomes))))
                (new-node       (add-node (copy-genome (first genomes))))
                (new-connection (add-connection (copy-genome (first genomes))))
                (copy           (copy-genome (first genomes))))
              new-population)))
    ; select a representative (from previous generation) and create new species' generation
    (make-instance 'species
                   :population new-population
                   :representative (alexandria:random-elt (population species)))))

(defmethod fitness ((species species))
  "Return a fitness value for the species as a whole"
  (let ((len (length (population species)))
        (fitness-list (mapcar #'fitness (population species))))
    (labels ((adjusted-fitness (f) (/ f len)))
      (if (= len 0)
          0.0
          (sum (mapcar #'adjusted-fitness fitness-list))))))

(defmethod reproduce ((generation generation) &rest rest)
  "Create a new generation from the previous generation."
  (declare (ignore rest))
  (let* ((species-list (mapcar #'initialize (species generation)))
         (species-fitness-list (mapcar #'fitness (species generation)))
         (total-fitness (sum species-fitness-list))
         ; Species population size is fractional by species fitness
         (species-population (mapcar (lambda (f) (floor ; species can go extinct if pop < 1
                                                   (* (population-size generation)
                                                      (/ f total-fitness))))
                                     species-fitness-list))
         ; Create a complete population of genomes based on species' fitness values
         (population (loop for species in (species generation)
                           for pop-size in species-population
                           nconc (reproduce species pop-size))))
    ; Speciate the population
    (setf species-list (speciate species-list population))
    ; Cull extinct species
    (setf species-list (remove-if (lambda (s)
                                    (= 0 (length (population s))))
                                  species-list))
    ; Return final genome
    (make-instance 'genome
                   :species species-list
                   :population-size (sum species-population))))


(defmethod update-network ((genome genome))
  "If required (as tracked by innovation number), update the network structure of the given genome.
   Returns a vector of nodes sorted by node id, containing (list id type value (list (next-node-pos weight) ...))"
  (when (>= (last-update genome) *innovation*)
    (return-from update-network (network genome)))
  (let* ((nodes (sort (nodes genome) #'< :key #'node-id))
         (node-count (length nodes))
         (network (make-array (list node-count) :initial-element nil)))
    ; Initialize network
    (loop for i from 0
          for node in nodes
          do (setf (aref network i)
                   (list (node-id node)
                         (node-type node)
                         0.0 ; value
                         nil))) ; initial list of connections
    ; Set connections for each node - enabled only.
    (dolist (conn (connections genome))
      (when (enabled conn)
        (let ((in-node (position (node-id (in-node conn)) 
                                 network
                                 :key #'first))
              (out-node (position (node-id (out-node conn))
                                  network
                                  :key #'first)))
          (push (list out-node (weight conn))
                (fourth (aref network in-node))))))
    ; Update the genome with this information
    (setf (last-update genome) *innovation*
          (network genome) network)))


(defmethod activate ((genome genome) inputs)
  "Activate the given genome and return the outputs.
   Activation continues for at most 100 cycles."
  ; Assign all nodes their respective set of connections
  ; Fire from input & bias nodes.
  ; Fire all subsequent nodes, until one of the following:
  ;     1. 100 such firings have occurred.
  ;     2. All output nodes have been touched.
  ; Return output nodes' current values.
  (let* ((network (update-network genome)) ; Gater all nodes that might be involved
         
         )
    
    )
  
  
  )
