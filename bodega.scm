(use numbers srfi-4 lapack-extras)

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
(define (z-p-g tf)
  (let ((num (car  tf))
	(den (cadr tf)))
    (values (poly-roots num) (poly-roots den) (/ (car num) (car den)))))

; evaluate the transfer function at each frequency
(define (freq-response tf freqs)
  (let-values (((zeros poles gain) (z-p-g tf)))
    (map
      (lambda (om)
	(let ((jw (make-rectangular 0 om)))
	  (* gain
	     (/ (foldl (lambda (acc x) (* acc (- jw x))) 1.0 zeros)
		(foldl (lambda (acc x) (* acc (- jw x))) 1.0 poles)))))
      freqs)))

; evaluate the bode plot for the transfer function
(define (bode tf start end step)
  ; convert x to dB
  (define ->dB
    (let ((ln10 (log 10)))
      (lambda (x) (* 20 (/ (log x) ln10)))))
  ; convert x to degrees
  (define ->deg
    (let ((k (/ 180 (atan 0 -1))))
      (lambda (x) (* k x))))

  #;(unless (proper? tf)
    (error "Improper transfer function"))

  (let* ((pulse (logspace start end step))
	 (resp  (freq-response tf pulse))
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
		; '((100 100) (1 110 1000))
		; '((1 0.1 7.5) (1 0.12 9 0 0))
		; '((-0.1 -2.4 -181 -1950) (1 3.3 990 2600))
		'((100) (1 30))
		-2 4 999)))
  (dump-gnuplot-script pulse mag phase))

; (print (z-p-g '((-10 20 0) (1 7 20 28 19 5))))
; (print (z-p-g '((100 100) (1 110 1000))))
; (print (poly-roots '(-10 20 0)))
