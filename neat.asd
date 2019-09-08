;;;; neat.asd

(asdf:defsystem #:neat
  :description "Describe neat here"
  :author "Jack Mudge <jakykong@theanythingbox.com>"
  :license  "LGPL v3+"
  :version "0.0.1"
  :serial t
  :depends-on (#:alexandria)
  :components ((:file "package")
               (:file "neat")))
