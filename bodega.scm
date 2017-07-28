(use numbers srfi-4 lapack-extras)

(define log10
  (let ((ln10 (log 10)))
    (lambda (x) (/ (log x) ln10))))

(define (proper? tf)
  (<= (length (car tf)) (length (cadr tf))))

(define (make-matrix m n)
  (make-f64vector (* m n) 0 #t))

(define (eval-poly poly x)
  (let loop ((deg (- (length poly) 1)) (coef poly) (acc 0.0))
    (if (zero? deg)
        (+ acc (car coef))
	(loop (- deg 1) (cdr coef) (+ acc (* (car coef) (expt x deg)))))))

; transform the polynomial into a monic one
(define (poly->monic-poly poly)
  (map (lambda (x) (/ x (car poly))) poly))

(define (companion poly)
  (let* ((ord (- (length poly) 1))
	 (mat (make-matrix ord ord)))
    ; zero-th grade terms in the latest column
    (do ((i (- ord 1) (+ i ord)) (c (reverse (cdr poly)) (cdr c)))
      ((>= i (* ord ord)))
      (f64vector-set! mat i (- (car c))))
    ; identity matrix shifted down by one row
    (do ((i ord (+ i (+ ord 1))))
      ((>= i (* ord ord)))
      (f64vector-set! mat i 1))
    mat))

(define dummy-vec (make-f64vector 0))

; find the polynomial roots by finding the eigenvalues of a matrix whose
; characteristic polynomial is 'poly'
(define (poly-roots* poly)
  (let ((mat (companion (poly->monic-poly poly)))
	(ord (- (length poly) 1)))
    (let-values
      (((a wr wi vl vr work)
	(dgeev "N" "N" ord mat
	       ; vr and vi
	       (make-f64vector ord 0.0)
	       (make-f64vector ord 0.0)
	       ; vl and vr aren't used
	       dummy-vec dummy-vec
	       ; workspace and its dimension
	       ; XXX: Evaluate the required size by passing -1 and fetch the
	       ; output in out[0]
	       (make-f64vector 16) 16)))
      (map make-rectangular (f64vector->list wr) (f64vector->list wi)))))

(define (poly-roots poly)
  ; do we have any term beside the constant term?
  (if (pair? (cdr poly))
      (poly-roots* poly)
      '()))

; sample n logarithmically-spaced values between 10^a and 10^b
(define (logspace a b n)
  (let ((step (/ (- b a) n)))
    (let loop ((n n) (e a) (v '()))
      (if (fx< n 0)
	  (reverse v)
	  (loop (fx- n 1) (+ e step) (cons (expt 10 e) v))))))

; transform the transfer function into its zero-pole-gain form
(define (zero-pole-gain tf)
  (let ((num (car  tf))
	(den (cadr tf)))
    (list (poly-roots num) (poly-roots den) (/ (car num) (car den)))))

; evaluate the transfer function at each frequency
(define (freq-response z+p+g freqs)
  (let-optionals z+p+g ((zeros 0) (poles 0) (gain 0))
    (map
      (lambda (om)
	(let ((jw (make-rectangular 0 om)))
	  (* gain
	     (/ (foldl (lambda (acc x) (* acc (- jw x))) 1.0 zeros)
		(foldl (lambda (acc x) (* acc (- jw x))) 1.0 poles)))))
      freqs)))

(define (dynamic z+p+g)
  (let* ((d+p   (map real-part (append (car z+p+g) (cadr z+p+g))))
	 (left  (- (apply max d+p)))
	 (right (- (apply min d+p))))
    ; cover the whole input range, if we have a single pole/zero then we center
    ; the whole plot around it
    ; XXX: pick an optimal number of divisions based on the [left,right] range
    (logspace (- (inexact->exact (round (log10 left)))  2)
	      (+ (inexact->exact (round (log10 right))) 2)
	      4000)))

; evaluate the bode plot for the transfer function
(define (bode tf #!optional pulse)
  ; convert x to dB
  (define (->dB x)
    (* 20 (log10 x)))
  ; convert x to degrees
  (define ->deg
    (let ((k (/ 180 (atan 0 -1))))
      (lambda (x) (* k x))))

  #;(unless (proper? tf)
    (error "Improper transfer function"))

  (let* ((z+p+g (zero-pole-gain tf))
	 (pulse (or pulse (dynamic z+p+g)))
	 (resp  (freq-response z+p+g pulse))
	 (phase (map (o ->deg angle) resp))
	 (mag   (map (o ->dB magnitude) resp)))
    (values pulse mag phase)))

; dump the data in a easy to plot format
(define (dump-gnuplot-script pulse mag phase)
  (with-output-to-file "mag.dat"
    (lambda ()
      (for-each
	(lambda (p v)
	  (printf "~a\t~a\n" (exact->inexact p) (exact->inexact v)))
	pulse mag)))
  (with-output-to-file "phase.dat"
    (lambda ()
      (for-each
	(lambda (p v)
	  (printf "~a\t~a\n" (exact->inexact p) (exact->inexact v)))
	pulse phase)))
  ; small script to setup the layout and plot the data
  (with-output-to-file "plot.gnuplot"
    (lambda ()
      (print
	(string-intersperse
	  '("unset key"
	    "set logscale x"
	    "set format x '10^{%L}'"
	    "set ytics 10"
	    "set multiplot layout 2, 1"
	    "plot 'mag.dat' with lines"
	    "set ytics 45"
	    "plot 'phase.dat' with lines"
	    "reset") "\n")))))

(let-values (((pulse mag phase)
	      (bode
		; system transfer function
		'((100 100) (1 110 1000))
		; '((1 0.1 7.5) (1 0.12 9 0 0))
		; '((-0.1 -2.4 -181 -1950) (1 3.3 990 2600))
		; '((100) (1 30))
		; frequency range
		#f #;(logspace -2 4 999)
		)))
  (dump-gnuplot-script pulse mag phase))

; (print (z-p-g '((-10 20 0) (1 7 20 28 19 5))))
; (print (z-p-g '((100 100) (1 110 1000))))
; (print (poly-roots '(-10 20 0)))
