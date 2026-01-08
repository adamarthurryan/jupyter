from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u

# from https://physics.mcmaster.ca/~harris/mwgc.dat

def load_clusters(): 

    gc_frame = Galactocentric()

    #load raw data from mwgc.dat.txt
    filename = "mwgc.dat.txt"
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

            clusters = []

            for row in lines[72:230]:
                if row.strip() :
                    name = row[0: 12].strip()
                    ra = row[26: 38].strip()
                    dec = row[39: 51].strip()
                    x = float(row[81:86].strip())-8.122  # adjust for sun position
                    y = float(row[86:92].strip())
                    z = float(row[92:].strip()) 

                    #print(x,y,z)
                    distance = float(row[69: 74].strip())
                    icrs = SkyCoord(ra=ra, dec=dec, distance=distance*u.kpc, unit=("hourangle", "deg"), frame='icrs')
                    galactocentric = SkyCoord(x=x*u.kpc, y=y*u.kpc, z=z*u.kpc, frame=gc_frame)
                    clusters.append({
                        "name": name,
                        "icrs": icrs,
                        "galactocentric": galactocentric,
                        "ra": ra,
                        "dec": dec,
                        "x": x,
                        "y": y,
                        "z": z,
                        "distance": float(distance)
                    })

            return clusters

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
