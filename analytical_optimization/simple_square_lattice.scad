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
 wb=2.62; va=2.69; lmax=3.60; hf=1.30; wa=0.60; vb=0.91; hc=0.20;description = "ana 3.6"; 
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
 wb=1.49; va=0.35; lmax=3.60; hf=0.50; wa=0.60; vb=3.25; hc=0.20; description = "ana broken 3.6";

w = wa+wb;
h = hf+hc;
l = va+vb;

tot_l = 50;
rep_z = 5;
rep_y = 5;

echo("interface size: ", w * rep_y, h * rep_z);
a = true; b = false;
              mat_a = a;

specimen_numbers = true;

brim = 10;

repeats = 5;
inter_dist = brim;

module finger()
{
    cube([l, wa, hf]);
}

module cell()
{
    finger();
    translate([vb,0,hf]) cube([va, w, hc]);
}
module pattern ()
{
    intersection()
    {
        cube([l,rep_y*w+wa, rep_z*h+hc]);
        union()
        {
            for (z = [-1:rep_z - 1])
                for (y = [0:rep_y])
                    translate([0,w*y,h*z+hc]) cell();
        }
    }
}

module sample(n=1)
{
    if (mat_a)
    {
        difference()
        {
            translate([-tot_l-brim,-brim,0]) cube([tot_l+brim+l,rep_y*w+wa+brim*2, .2]);
            translate([0,0,-1])cube([tot_l,rep_y*w+wa, rep_z*h+hc]);
        }
        pattern();
        difference()
        {
            translate([-tot_l,0,0]) cube([tot_l,rep_y*w+wa, rep_z*h+hc]);
            translate([-tot_l + rep_y*w / 2, rep_y*w / 2, rep_z*h+hc-.1]) linear_extrude(1.0) text(str(n),halign="center", valign="center",size = 7*10/(rep_y*w));
        }
    }
    else
    {
        translate([l,-brim,0]) cube([tot_l+brim,rep_y*w+wa+brim*2, .2]);
        difference()
        {
            cube([tot_l+l,rep_y*w+wa, rep_z*h+hc]);
            translate([tot_l - rep_y*w / 2, rep_y*w / 2, rep_z*h+hc-.1]) linear_extrude(1.0) text(str(n),halign="center", valign="center",size = 7*10/(rep_y*w));
            pattern();
        }
    }
}

for (y = [0:repeats-1])
{
    translate([0, y * (rep_y*w+wa + inter_dist) ,0]) sample(y + 1);
}
if (mat_a)
{
    difference()
    {
        translate([-tot_l,-12-brim,0]) cube([tot_l,12,.3]);
        translate([-4,-6-brim,-.1]) linear_extrude([.5])
        rotate([0,0,180]) text(description, valign="center",size=5);
    }
}
else
{
    difference()
    {
        translate([l+0,-12-brim,0]) cube([tot_l,12,.3]);
        translate([l+4,-6-brim,-.1]) linear_extrude([.5])
        text(description, valign="center",size=5);
    }
}