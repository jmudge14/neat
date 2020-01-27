;;;; neat.lisp

(in-package #:neat)

(import 'anaphora:aif)
(import 'anaphora:it)

;;;; Parameters
(defparameter +connection-perturb-probability+ 0.90
  "Chance of any particular connection being perturbed during perturbance mutation" )
(defparameter +connection-perturb-scale+ 0.40
  "Scale factor applied to a random number in (-1.0, 1.0) to perturb weights by")

(defparameter +new-connection-weight-magnitude+ 1.0
  "Largest absolute value of new connections' weights, selected randomly between (-mag,+mag)")

; Coeifficients contributing to the distance calculation.
(defparameter +distance-excess-coefficient+ 1.0
  "Coefficient times number of excess genes")
(defparameter +distance-disjoint-coefficient+ 1.0
  "Coefficient times number of disjoint genes")
(defparameter +distance-weights-coefficient+ 2.0
  "Coefficient times average difference in weights")

(defparameter +compatibility-distance-threshold+ 1.0 
  "Maximum distance for two genomes to be of the same species.")


(defparameter +equal-fitness-threshold+ 0.2
  "Difference in fitness for parents to be considered equal during breeding")

(defparameter +default-population-size+ 150
  "The default population size of generations in a new experiment.")

(defparameter +species-mutation-rates+
  '((breed          . 2000)  ; (breed g1 g2)
    (weight         . 1000)  ; (perturb g1)
    (new-node       .   70)  ; (add-node g1)
    (new-connection .  200)  ; (add-connection g1)
    (copy           .    0)) ; (copy-genome g1) ; Set to 0 because carry-over percent supersedes it.
  "Relative weight of various mutation types. 
   Copy does not mutate. Breed is sexual reproduction.
   All others are asexual reproduction.")

(defparameter +generation-carryover-percent+ 0.10
  "Percentile of best performers carried over from previous generation")

(defparameter +max-activate-iterations+ 10
  "Maximum number of firing rounds during activation if all output nodes are not yet reached.
   If this value is larger, deeper networks will fully activate, but activation could take longer.")



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
    (push node *nodes*)
    node))


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

(defmethod perturb ((genome genome))
  "Randomly adjust connection weights"
  (dolist (connection (connections genome) genome)
    (let* ((new-weight (+ (weight connection)
                          (* (alexandria:gaussian-random -1.0 1.0)
                             +connection-perturb-scale+)))
           (perturb? (> +connection-perturb-probability+ 
                        (random 1.0))))
      (when perturb?
        (setf (weight connection) new-weight)))))

#| (defmacro where (&rest l) `(lambda (it) (or ,@l))) |#

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
   (population-size :initform +default-population-size+ :initarg :population-size :accessor population-size)))


(defun make-generation (genome-list &optional (species-list nil) (population-size +default-population-size+))
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
    #|
    (make-instance 'species
                   :population new-population
                   :representative (alexandria:random-elt (population species))) |#
    new-population
    ))

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
                                                      (if (plusp total-fitness) 
                                                          (/ f total-fitness)
                                                          0.0))))
                                     species-fitness-list))
         ; Create a complete population of genomes based on species' fitness values
         (population (loop for species in (species generation)
                           for pop-size in species-population
                           nconc (reproduce species pop-size))))
    ; Include the most fit 10% of individuals from the previous generation
    ; This causes a slight over-size on each generation compared to target, but ensures good 
    ; performers stay in the gene pool.
    (let* ((num-to-keep (floor (* +generation-carryover-percent+ (population-size generation))))
           (prev-pop (sort (population generation) #'> :key #'fitness)))
      (setf population
            (nconc population 
                   (if (> num-to-keep (length prev-pop))
                       prev-pop
                       (subseq prev-pop 0 num-to-keep)))))
    ; Speciate the population
    (setf species-list (speciate species-list population))
    ; Cull extinct species
    (setf species-list (remove-if (lambda (s)
                                    (= 0 (length (population s))))
                                  species-list))
    ; Return final genome
    (make-instance 'generation
                   :species species-list
                   ; :population-size (sum species-population)
                   :population-size (population-size generation)
                   )))


(defmethod update-network ((genome genome))
  "If required (as tracked by innovation number), update the network structure of the given genome.
   Returns a vector of nodes sorted by node id, containing (list id type value (list (next-node-pos weight) ...) touched)"
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
                         nil ; initial list of connections 
                         nil ; touched flag - only useful for output nodes.  
                         ))) 
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

(defun sigmoid (num)
  (/ 1 (+ 1 (exp (- num)))))


(defmethod activate ((genome genome) inputs)
  "Activate the given genome and return the outputs.
   Activation continues for at most 100 cycles."
  ; Assign all nodes their respective set of connections
  ; Fire from input & bias nodes.
  ; Fire all subsequent nodes, until one of the following:
  ;     1. +max-activate-iterations+ such firings have occurred.
  ;     2. All output nodes have been touched.
  ; Return output nodes' current values.
  (let* ((network (update-network genome)) ; Gather all nodes that might be involved
         (out-nodes (loop for i from 0 upto (1- (array-dimension network 0)) ; List of output node indices
                          when (eq (second (aref network i)) :output) collect i)))
    ; Reset initial values (note - assumes feed-forward expectation
    (loop for i from 0 upto (1- (array-dimension network 0)) do
          (setf (third (aref network i))
                (case (second (aref network i))
                  (:sensor (pop inputs))
                  (:bias 1.0)
                  (t 0.0))))
    (block firing-loop 
       (dotimes (iter +max-activate-iterations+)
         ; Fire all nodes.
         (loop for node being the elements of network do
               (let ((cur-value (tanh (third node)))
                     (connections (fourth node)))
                 (loop for conn in connections do
                       (destructuring-bind (out-node-id weight) conn 
                         (let ((out-node (aref network out-node-id)))
                           (incf (third out-node)
                                 (* cur-value weight)) ; TODO use transfer function?
                           ; Update touched if it is an output node
                           (when (eq (second out-node) :output)
                             (setf (fifth out-node) t))))))
               ; Reset current node value for sensor and hidden nodes, 
               (unless (find (second node) '(:output)) 
                 (setf (third node) 0.0)))
         ; Check if we have completed firing before +max-activate-iterations+
         ; That is, if all output nodes have been reached.
         (loop for i in out-nodes
               with all = t
               unless (fourth (aref network i)) 
               do (setf all nil) ; mark not all if anything not reached
               end
               finally (when all (return-from firing-loop)))))
    ; Return the final output results
    (loop for i in out-nodes
          collect (third (aref network i)))))

;;; TODO

(defun make-genome (sensor-count output-count)
  "Return a fully connected, randomized, initial genome with the given number of sensors and outputs.
   This is primarily needed to create the first generation of an experiment."
  (let ((genome (make-instance 'genome))
        (in-nodes nil)
        (out-nodes nil))
    ; Create nodes
    (add-node-of-type genome :bias)
    (dotimes (n sensor-count)
      (add-node-of-type genome :sensor))
    (dotimes (n output-count)
      (add-node-of-type genome :output))
    ; Populate lists of nodes to be interconnected
    (loop for node in (nodes genome) do
          (if (eq (node-type node) :output)
              (push node out-nodes)
              (push node in-nodes)))
    ; Connect every sensor node to every output node
    (loop for in-node in in-nodes do
          (loop for out-node in out-nodes do
                (add-connection genome (list in-node out-node))))
    ; Return genome for further consumption.
    genome))

(defun make-first-generation (sensor-count output-count &optional (population-size +default-population-size+))
  "Return the first generation of a new experiment."
  (let* ((starting-genome (make-genome sensor-count output-count))
         (genome-list (list starting-genome)))
    ; Create initial population of genomes with random mutations
    (dotimes (n (1- population-size))
      (let ((next-genome (copy-genome starting-genome)))
        (perturb next-genome)
        (push next-genome genome-list)))
    ; Return a generation from the above initial population
    (make-generation genome-list nil population-size)))


(defmethod population ((generation generation))
  "Return all genomes present in a generation."
  (loop for species in (species generation)
        appending (population species)))

(defmethod max-fitness ((generation generation))
  (loop for genome in (population generation)
        maximizing (fitness genome)))

;;;; XOR experiment

(defparameter +expected-values-xor+ 
  '(((0 0) 0.0)
    ((0 1) 1.0)
    ((1 0) 1.0)
    ((1 1) 0.0)))

(defmethod xor-trial-fitness ((genome genome))
  "Tests a genome against a set of expected inputs and outputs.
   Assumes a single output is required, for simplicity."
  (loop for trial in +expected-values-xor+
        for inputs = (first trial) then (first trial)
        for target = (second trial) then (second trial)
        for outval = (first (activate genome inputs))
        with fitness = 0.0
        do (incf fitness (abs (- outval target)))
        finally (return (- 10.0 fitness))))

(defmethod update-fitness-singles ((generation generation) fitness-func)
  (loop for genome in (population generation)
        do (setf (fitness genome) (funcall fitness-func genome))))


(defmethod get-next-generation-single-fitness ((generation generation) fitness-func)
  (update-fitness-singles generation fitness-func)
  (let ((next-gen (reproduce generation)))
    (update-fitness-singles next-gen fitness-func)
    ; Return the previous generation if we didn't improve overall fitness.
    (values next-gen (max-fitness next-gen))))


(defun run-xor-experiment ()
  (let ((gen (make-first-generation 2 1 100)))
    (dotimes (n 1000) 
      (setf gen (get-next-generation-single-fitness gen #'xor-trial-fitness))
      (format t "Trial ~A - Max ~A~%" n (max-fitness gen))
      (when (> (max-fitness gen) 9.9)
        (format t "Took ~A trials to reach fitness 9.9." n)
        (return-from run-xor-experiment))))
  (format t "Could not reach target within 1000 trials."))
