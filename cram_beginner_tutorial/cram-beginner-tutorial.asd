(defsystem cram-beginner-tutorial
  :depends-on (roslisp cram-language turtlesim-msg cl-transforms geometry_msgs-msg designators cram-reasoning 
                 cram-language-designator-support actionlib actionlib_tutorials-msg process-modules turtle_actionlib-msg cl-tf)
  :components
  ((:module "src"
            :components
            ((:file "package")
             (:file "control-turtlesim" :depends-on  ("package"))
             (:file "simple-plans" :depends-on  ("package" "control-turtlesim"))
             (:file "action-designators" :depends-on  ("package"))
             (:file "turtle-action-client" :depends-on  ("package"))
             (:file "location-designators" :depends-on  ("package"))
             (:file "process-modules" :depends-on  ("package" "control-turtlesim" "simple-plans" "action-designators" "turtle-action-client" "location-designators"))))))

