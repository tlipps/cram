;;;
;;; Copyright (c) 2018, Christopher Pollok <cpollok@uni-bremen.de>
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are met:
;;;
;;;     * Redistributions of source code must retain the above copyright
;;;       notice, this list of conditions and the following disclaimer.
;;;     * Redistributions in binary form must reproduce the above copyright
;;;       notice, this list of conditions and the following disclaimer in the
;;;       documentation and/or other materials provided with the distribution.
;;;     * Neither the name of the Institute for Artificial Intelligence/
;;;       Universitaet Bremen nor the names of its contributors may be used to
;;;       endorse or promote products derived from this software without
;;;       specific prior written permission.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
;;; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
;;; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
;;; LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
;;; CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
;;; SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
;;; CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
;;; ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;; POSSIBILITY OF SUCH DAMAGE.

(in-package :pr2-em)

;;; OBJECT-INTERFACE METHODS

(defparameter *drawer-handle-grasp-x-offset* -0.02 "in meters")
(defparameter *drawer-handle-pregrasp-x-offset* 0.10 "in meters")
(defparameter *drawer-handle-lift-x-offset* 0.48 "in meters")
(defparameter *drawer-handle-2nd-lift-x-offset* 0.58 "in meters")

; Might be necessary to find out what kind of handle we are dealing with. But we could also just open wide and be done with it.
(defmethod obj-int:get-object-type-gripper-opening ((object-type (eql :container)))
  0.09)

; Find out where the handle is and calculate the transform from there.
(defmethod obj-int:get-object-type-to-gripper-transform ((object-type (eql :container))
                                                         object-name
                                                         arm
                                                         (grasp (eql :front)))
  (setf object-name (roslisp-utilities:rosify-underscores-lisp-name object-name))
  (let* ((handle-name (cl-urdf:name (get-handle-link object-name)))
         (handle-tf (cl-tf:transform->transform-stamped "map" handle-name 0
                                                        (cl-tf:pose->transform
                                                         (get-urdf-link-pose handle-name))))
         (container-tf (cl-tf:transform->transform-stamped "map" object-name 0
                                                           (cl-tf:pose->transform
                                                            (get-urdf-link-pose object-name))))
         (tool-frame (ecase arm
                       (:left cram-tf:*robot-left-tool-frame*)
                       (:right cram-tf:*robot-right-tool-frame*))))
     (cram-tf:multiply-transform-stampeds object-name
                                          tool-frame
                                          (cram-tf:multiply-transform-stampeds object-name handle-name
                                                                               (cram-tf:transform-stamped-inv container-tf)
                                                                               handle-tf)
                                          (cl-transforms-stamped:make-transform-stamped
                                           handle-name
                                           tool-frame
                                           0.0
                                           (cl-transforms:make-3d-vector *drawer-handle-grasp-x-offset* 0.0d0 0.0d0)
                                           (cl-transforms:matrix->quaternion
                                            #2A((0 0 -1)
                                                (0 1 0)
                                                (1 0 0)))))))

; Should be fine without a joint-type.
(defmethod obj-int:get-object-type-to-gripper-pregrasp-transform ((object-type (eql :container))
                                                                   object-name
                                                                   arm
                                                                   (grasp (eql :front))
                                                                   grasp-transform)
  (cram-tf:translate-transform-stamped grasp-transform :x-offset *drawer-handle-pregrasp-x-offset*))

; We need the joint-type, maybe add a optional parameter to the generic.
; But for now I can live with this.
(defmethod obj-int:get-object-type-to-gripper-lift-transform ((object-type (eql :container))
                                                              object-name
                                                              arm
                                                              (grasp (eql :front))
                                                              grasp-pose)
  (let ((grasp-pose grasp-pose))
    (cram-tf:translate-transform-stamped grasp-pose :x-offset *drawer-handle-lift-x-offset*)))

(defmethod obj-int:get-object-type-to-gripper-2nd-lift-transform ((object-type (eql :container))
                                                                  object-name
                                                                  arm
                                                                  (grasp (eql :front))
                                                                  grasp-pose)
  (cram-tf:translate-transform-stamped grasp-pose :x-offset *drawer-handle-2nd-lift-x-offset*))
