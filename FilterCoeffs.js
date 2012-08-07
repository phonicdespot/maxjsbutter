//================================================================
//
// Calculate 4th order Butterworth low pass filter coefficients
//
//================================================================

inlets = 4;     /*  This gives the js object two inlets.

                    The first inlet is used to change the
                    cut-off frequency.
                    
                    The second inlet is used to change the
                    sampling frequency.

		    The third inlet sets the filter order.

		    The fourth inlet sets the filter type
		    (low- or high-pass)
                */

setinletassist(0,"Cut-off frequency");
setinletassist(1,"Sampling frequency");
setinletassist(2,"Filter order");
setinletassist(3,"Filter type");


outlets = 1;    /*  This gives the js object a single outlet.

                    This is used to output the list to set the
                    coefficients of a Max/MSP biquad~ object.
                */


//________________________________________________________________
//
// Global variables that can be changed using the inlets to the js
// object
//                
var fs = 44100;             // The sampling frequency (Hz)
var fc = 440;               // The cut-off frequency (Hz)


//________________________________________________________________
//
// This is the prototype 3rd order Butterworth LPF (i.e. the
// cut-off angular frequency is 1).
//
// The num field is an array containing the coefficients of the
// numerator of the transfer function.
//
// The den field is an array containing the coefficients of the
// denominator of the transfer function.
//
// This is constant.
//
var butter2 = {
    num: [0,0,1],
    den: [1,Math.sqrt(2),1]
};

var butter3 = new Array(2);
butter3[0] = {
    num: [0,0,1],
    den: [1,1,1]
};
butter3[1] = {
    num: [0,0,1],
    den: [0,1,1]
};

var butter4 = new Array(2);
butter4[0] = {
    num: [0,0,1],
    den: [1,2*Math.cos(Math.PI/8),1]
};
butter4[1] = {
    num: [0,0,1],
    den: [1,2*Math.cos(3*Math.PI/8),1]
};


//________________________________________________________________
//
// These are the remaining global variables 
//
var T = 1/fs;               // The sampling period
var wc = 2*Math.PI*fc;      // The cut-off angular frequency (rad/s)
var pwc = prewarp(wc);      // The prewarped cut-off angular freq.

var lpf_coeffs = lpf(butter4,pwc);  // The coefficients of the
                                    // analogue filter that is to be
                                    // transformed using the bilinear
                                    // transformation.

var coeffs = blt(lpf_coeffs,T);     // The coefficients of the digital
                                    // filter.




//================================================================
//
// Public functions
//
//================================================================


//________________________________________________________________
//
//  This responds to integers in the inlets
//
function msg_int(a) {
	switch (inlet) {
		case 0 :
		    // Change the cut-off frequency
		    fc = a;
			break;
		case 1 :
		    // Change the sampling frequency
			fs = a;
	}
	calculate_coeffs();
    output_coeffs();
}

//________________________________________________________________
//
//  This responds to floats in the inlets
//
function msg_float(a) {
	switch (inlet) {
		case 0 :
		    // Change the cut-off frequency
		    fc = a;
	}
	calculate_coeffs();
    output_coeffs();
}

//________________________________________________________________
//
//  Output the biquad coefficients when a bang is received.
//
function bang() {
    output_coeffs();
}

//________________________________________________________________
//
//  Output the biquad coefficients on loading.
//
function loadbang() {
    output_coeffs();
}


//================================================================
//
// Private functions
//
//================================================================


//________________________________________________________________
//
//  Output the biquad coefficients.
//
output_coeffs.local = 1;    // Make the function local (private)
function output_coeffs() {
    var bq_cf = convert_to_cascade_list(coeffs);
    outlet(0,bq_cf);
    
    // Post the coefficients in the Max window
    post('cascade coefficients',bq_cf,'\n');
}

//________________________________________________________________
//
//  Calculate the coefficients when either the cut-off frequency
//  or the sampling frequency has changed.
//
calculate_coeffs.local = 1;
function calculate_coeffs() {
    T = 1/fs;
    wc = 2*Math.PI*fc;
    pwc = prewarp(wc);
    lpf_coeffs = lpf(butter4,pwc);
    coeffs = blt(lpf_coeffs,T);
}

//________________________________________________________________
//
//  Prewarp the cut-off (angular) frequency when designing the
//  analogue filter before applying the bilinear transform.
//
prewarp.local = 1;
function prewarp(w) {
    return (2/T)*Math.atan(w*T/2);
}

//________________________________________________________________
//
//  Concatanate the lists of biquad coefficients.
//
function convert_to_cascade_list(cf) {
    var a = [];
    for (var i=0; i<2; i++) {
        a.push.apply(a,convert_to_biquad_list(cf[i]));
    }
    return a;
}

//________________________________________________________________
//
//  Convert the coefficients from the num/den form used in
//  butter2, etc. to a list that can be used by the Max/MSP
//  biquad~ object.
//
//  The leading coefficient in the denominator must be 1,
//  so the coefficients are scaled accordingly.
//
//  I wondered about using a loop here but decided it wasn't
//  worth it.
//
convert_to_biquad_list.local = 1;
function convert_to_biquad_list(cf) {
    var result = new Array(5);
    result[0] = cf.num[0]/cf.den[0];
    result[1] = cf.num[1]/cf.den[0];
    result[2] = cf.num[2]/cf.den[0];
    result[3] = cf.den[1]/cf.den[0];
    result[4] = cf.den[2]/cf.den[0];
    return result;
}

//________________________________________________________________
//
//  This function takes the coefficients for the prototype analogue
//  low pass filter and transforms them into the coefficients for
//  an analogue low pass filter with the specified cut-off frequency.
//
lpf.local = 1;
function lpf(cf,w) {
    var result = new Array(2);
    for (var i=0; i<2; i++) {
        result[i] = {
            num: lpf_qd(cf[i].num,w),
            den: lpf_qd(cf[i].den,w)
        }
    }
    return result;
}

//________________________________________________________________
//
//  Given a 3 element array [a,b,c] and a scalar w, this function
//  returns the array [a,b*w,c*w^2].
//
//  This is only used in the lpf function.
//
lpf_qd.local = 1;
function lpf_qd(cf,w) {
    var result = new Array(3);
    for (var i=0; i<3; i++) {
        result[i] = cf[i]*Math.pow(w,i);
    }
    return result;
}

//________________________________________________________________
//
//  This function takes the coefficients for an analogue filter
//  and uses the bilinear transform to generate an equivalent
//  digital filter.
//
blt.local = 1;
function blt(cf,T) {
    var result = new Array(2);
    for (var i=0; i<2; i++) {
        result[i] = {
            num: blt_qd(cf[i].num,T),
            den: blt_qd(cf[i].den,T)
        }
    }
    return result;
}
blt_qd.local = 1;
function blt_qd(cf,T) {
    var result = new Array(3);
    result[0] = blt_qd_cf(cf,[4,2,1],T);
    result[1] = blt_qd_cf(cf,[-8,0,2],T);
    result[2] = blt_qd_cf(cf,[4,-2,1],T);
    return result;
}
blt_qd_cf.local = 1;
function blt_qd_cf(cf,cf2,T) {
    var result = 0;
    for (var i=0; i<3; i++) {
        result += cf[i]*cf2[i]*Math.pow(T,i);
    }
    return result;
}

