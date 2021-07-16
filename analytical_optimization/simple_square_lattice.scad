////fem:
// wb=2.70; va=2.88; lmax=3.60; hf=0.80; wa=0.60; vb=0.72; hc=0.20; F=22.56;
// wb=2.10; va=1.68; lmax=2.40; hf=0.80; wa=0.60; vb=0.72; hc=0.20;
////from fem draft2:
// wb=2.54; va=2.82; lmax=3.6; hf=0.8; wa=0.6; vb=lmax-va; hc=0.2;
// wb=2.35; va=2.27; lmax=3.0; hf=0.8; wa=0.6; vb=lmax-va; hc=0.2;
// wb=2.21; va=1.84; lmax=2.4; hf=0.7; wa=0.6; vb=lmax-va; hc=0.2;
// wb=2.01; va=1.36; lmax=1.8; hf=0.6; wa=0.6; vb=lmax-va; hc=0.2;
//// ana correct:
// old print wb=2.45; va=2.78; lmax=3.60; hf=1.00; wa=0.60; vb=0.82; hc=0.20;
// wb=2.68; va=2.66; lmax=3.60; hf=1.09; wa=0.60; vb=0.94; hc=0.20;
// wb=2.69; va=2.22; lmax=3.00; hf=0.91; wa=0.60; vb=0.78; hc=0.20;
// wb=2.67; va=1.78; lmax=2.40; hf=0.73; wa=0.60; vb=0.62; hc=0.20;
// wb=2.30; va=1.27; lmax=1.80; hf=0.62; wa=0.60; vb=0.53; hc=0.20;
// printable:
// wb=2.62; va=2.69; lmax=3.60; hf=1.30; wa=0.60; vb=0.91; hc=0.20;description = "ana 3.6"; 
// wb=2.57; va=2.23; lmax=3.00; hf=1.10; wa=0.60; vb=0.77; hc=0.20;description = "ana 3.0"; 
// wb=2.51; va=1.78; lmax=2.40; hf=0.90; wa=0.60; vb=0.62; hc=0.20; description = "ana 2.4"; 
// wb=2.41; va=1.32; lmax=1.80; hf=0.70; wa=0.60; vb=0.48; hc=0.20; description = "ana 1.8"; 
//// ana correct broken:
// old print wb=2.45; va=1.15; lmax=3.60; hf=1.00; wa=0.60; vb=2.45; hc=0.20;
//wb=1.52; va=0.34; lmax=3.60; hf=0.49; wa=0.60; vb=3.26; hc=0.20;
// wb=1.38; va=0.29; lmax=3.00; hf=0.45; wa=0.60; vb=2.71; hc=0.20;
// wb=1.22; va=0.30; lmax=2.40; hf=0.39; wa=0.60; vb=2.10; hc=0.20;
// wb=1.04; va=0.30; lmax=1.80; hf=0.33; wa=0.60; vb=1.50; hc=0.20;
// printable:
// wb=1.49; va=0.35; lmax=3.60; hf=0.50; wa=0.60; vb=3.25; hc=0.20; description = "ana broken 3.6"; 

lmax = 3.6;

sample_points =
    [
    [[2.4, 2.7, 3.6, 0.8], "a"], // >> near FEM opt
    [[2.1, 2.7, 3.6, 0.8], "b"], // wb-
    [[2.7, 2.7, 3.6, 0.8], "c"], // wb+
    [[2.4, 2.4, 3.6, 0.8], "d"], // va-  *
    [[2.4, 3.0, 3.6, 0.8], "e"], // va+
    [[2.4, 2.7, 3.6, 0.6], "f"], // hf-  *
    [[2.4, 2.7, 3.6, 1.0], "g"], // hf+
    [[1.5, 0.3, 3.6, 0.5], "h"], // broken
    [[1.2, 0.3, 3.6, 0.5], "i"], // wb-
    [[1.7, 0.3, 3.6, 0.5], "j"], // wb+
    [[1.5, 0.6, 3.6, 0.5], "k"], // va+
    [[1.5, 0.3, 3.6, 0.3], "l"], // hf-
    [[1.5, 0.3, 3.6, 0.7], "m"], // hf+
//    [[2.1, 2.4, 3.6, 0.8], ""], // wb- from d
//    [[2.7, 2.4, 3.6, 0.8], ""], // wb+ from d
//    [[2.4, 2.1, 3.6, 0.8], ""], // va- from d
//  //[[2.4, 2.7, 3.6, 0.8], "b"], // va+ from d
//    [[2.4, 2.4, 3.6, 0.6], ""], // hf- from d
//    [[2.4, 2.4, 3.6, 1.0], ""], // hf+ from d
//    [[2.7, 2.7, 3.6, 0.6], ""], // wb- from f
//    [[2.1, 2.7, 3.6, 0.6], ""], // wb+ from f
//    [[2.4, 2.4, 3.6, 0.6], ""], // va- from f
//    [[2.4, 3.0, 3.6, 0.6], ""], // va+ from f
//    [[2.4, 2.7, 3.6, 0.4], ""], // hf- from f
//  //[[2.4, 2.7, 3.6, 1.0], "g"], // hf+ from f
//    [[2.4, 2.4, 3.6, 0.6], ""], // near FEM opt va- & hf-
//    [[2.1, 2.4, 3.6, 0.6], ""], // wb-
//    [[2.7, 2.4, 3.6, 0.6], ""], // wb+
//    [[2.4, 2.1, 3.6, 0.6], ""], // va-
//  //[[2.4, 2.7, 3.6, 0.6], "f"], // va+
//    [[2.4, 2.4, 3.6, 0.4], ""], // hf-
//  //[[2.4, 2.4, 3.6, 0.8], "d"], // hf+
    ];

tpla = true; pp = false; mat_a = tpla;

n=2;

N = len(sample_points);

use_brim = false;
brim = 7;
outer_brim = 20;
side_brim = 30;

tot_l = 50;
rep_z = 5;
rep_y = 5;

specimen_numbers = true;


wmin = .3;
hmin = .2;

repeats = 5;

module finger(wa,wb,va,vb,hc,hf,w,h,l)
{
    cube([l, wa, hf]);
}

module cell(wa,wb,va,vb,hc,hf,w,h,l)
{
    finger(wa,wb,va,vb,hc,hf,w,h,l);
    translate([vb,0,hf]) cube([va, w, hc]);
}
module pattern (wa,wb,va,vb,hc,hf,w,h,l)
{
    intersection()
    {
        cube([l,rep_y*w+wa, rep_z*h+hc]);
        union()
        {
            for (z = [-1:rep_z - 1])
                for (y = [0:rep_y])
                    translate([0,w*y,h*z+hc]) cell(wa,wb,va,vb,hc,hf,w,h,l);
        }
    }
}


module sample(tag, wb,va,lmax,hf, extend_brim_after)
{
    sample_(tag,
        2 * wmin,wb,
        va,lmax - va,
        hmin,hf,
        wmin+wb,hmin+hf,lmax,
        extend_brim_after);
}
module sample_(tag, wa,wb,va,vb,hc,hf,w,h,l, extend_brim_after)
{
    * echo("interface size: ", w * rep_y, h * rep_z);

    if (mat_a)
    {
        if (use_brim) 
        { // brim
            difference() 
            {
                translate([-tot_l-side_brim,-brim,0])
                    cube([tot_l+side_brim+l,rep_y*w+wa+brim*2, .2]);
                translate([0,0,-1])
                    cube([tot_l,rep_y*w+wa, rep_z*h+hc]);
            }
            if (extend_brim_after)
                translate([-tot_l-side_brim,rep_y*w,0])
                cube([tot_l+side_brim+l,outer_brim, .2]);
        }
        pattern(wa,wb,va,vb,hc,hf,w,h,l);
        difference()
        {
            translate([-tot_l,0,0]) cube([tot_l,rep_y*w+wa, rep_z*h+hc]);
            translate([-tot_l + 3, rep_y*w / 2, rep_z*h+hc-.1]) linear_extrude(1.0) text(tag,halign="left", valign="center",size = min(5,6*.1*(rep_y*w)));
        }
    }
    else
    {
        if (use_brim) 
        { // brim
            translate([l,-brim,0])
            cube([tot_l+side_brim,rep_y*w+wa+brim*2, .2]);
            if (extend_brim_after)
                translate([l,rep_y*w,0])
                cube([tot_l+side_brim,outer_brim, .2]);
        }
        difference()
        {
            cube([tot_l+l,rep_y*w+wa, rep_z*h+hc-.0001]);
            translate([tot_l - 0, rep_y*w / 2, rep_z*h+hc-.1]) linear_extrude(1.0) text(tag,halign="right", valign="center",size = min(5,6*.1*(rep_y*w)));
            pattern(wa,wb,va,vb,hc,hf,w,h,l);
        }
    }
}


module samples(from, till)
{
    sample(str(sample_points[(from + N - round((n-1) / 5 * N)) % N][1], n)
        , sample_points[(from + N - round((n-1) / 5 * N)) % N][0][0]
        , sample_points[(from + N - round((n-1) / 5 * N)) % N][0][1]
        , sample_points[(from + N - round((n-1) / 5 * N)) % N][0][2]
        , sample_points[(from + N - round((n-1) / 5 * N)) % N][0][3], from == till);
    if (till > from)
    {
        translate([0, brim + rep_y*(sample_points[from][0][0]+2*wmin),0])
        samples(from + 1, till);
    }
}

*samples(0, 1);
samples(0, len(sample_points) - 1);

if (mat_a)
{
    translate([-tot_l-side_brim,-outer_brim,0])
    if (use_brim) cube([tot_l+side_brim+lmax,outer_brim, .2]);
}
else
{
    translate([lmax,-outer_brim,0])
    if (use_brim) cube([tot_l+side_brim,outer_brim, .2]);
}
