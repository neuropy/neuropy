''' colortable used for displaying rf array obtained by reverse correlation'''

from Core import *

def rcLUT(n=256):
    ''' n is the number of entries you want in the returned colortable. Defaults to 256.'''

    v = array(range(n-1,-1,-1),dtype=float,ndmin=2).transpose()/max(n-1,1) # normal colormap. Red is high value, blue is low, green is medium
    #v = array(range(0,n,1),ndim=2).transpose()/max(n-1,1); # inverted colormap
    m = 64.0/255.0

    # fast way

    r = np.exp(-(v/m)**2 / 2)
    g = np.exp(-((v-127.0/255.0)/m)**2 / 2)
    b = np.exp(-((v-1)/m)**2 / 2)

    colortable = cat((r, g, b),1)
    return colortable

    # slow way

    #colortable = []
    #for i in range(1,n+1,1)
    #
    #	r = v(i);
    #	g = v(i)-127/255;
    #	b = v(i)-1;
    #	colortable = [colortable; exp(-(r/m)^2 / 2)  exp(-(g/m)^2 / 2) exp(-(b/m)^2 / 2) ];
    #
    #end;

    # Delphi code

    #function TMSeqForm.Grey2Color (gval : smallint) : TColor;
    #var gr,gg,gb,g : smallint;
    #begin
    #  g := 255-gval;
    #  gr:= g;
        #  gg:= g-127;
    #  gb:= g-255;
    #
    #  Result:= RGB(round(255*exp(-sqr(gr/64)/2)),  {Color}
    #               round(255*exp(-sqr(gg/64)/2)),
    #               round(255*exp(-sqr(gb/64)/2)));
