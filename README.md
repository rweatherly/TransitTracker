# TransitTracker

Finds ISS and Tiangong transits of the Sun and Moon for one observing site.

`transit_finder.py` scans the next 14 days and reports:

- **TRANSIT**: the station's path crosses the disc as seen from the site.
- **POSSIBLE**: a near miss (within one disc radius of the limb) three or more
  days out, where orbit-prediction drift could still move the path over the site.

## Observing site

The site is never committed. Give it at run time, in any of these ways:

```
python transit_finder.py --lat 40.0 --lon -80.0 --elev 180
TT_LAT=40.0 TT_LON=-80.0 TT_ELEV=180 python transit_finder.py
cp location.example.json location.json   # edit; location.json is gitignored
```

## Orbit data

The GitHub Action in `.github/workflows/fetch-tle.yml` pulls fresh ISS and
Tiangong TLEs from CelesTrak once a day into `data/stations.tle`. The script
uses that file when present, otherwise it downloads from CelesTrak directly.
`--tle FILE` overrides both. It warns when the TLEs are over three days old.

## Setup

```
pip install -r requirements.txt
```

`skyfield-data` bundles the de421 ephemeris. Without it, Skyfield downloads
de421 into `cache/` on first run.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--days` | 14 | how far ahead to scan |
| `--min-alt` | 10 | skip events with the Sun/Moon below this altitude |
| `--low-alt` | 20 | mark events below this altitude as `low` (may be blocked) |
| `--bodies` | sun moon | which discs to check |
| `--possible-after-days` | 3 | flag near misses only this far out or more |
| `--possible-margin` | 1.0 | near-miss window, in disc radii past the limb |
| `--tz` | America/Detroit | time zone for local times |
| `--start` | now | scan start, ISO UTC |
| `--json` | off | machine-readable output |

## Output fields

`direction` is the compass point of the event and `low` is true when the Sun or
Moon is under `--low-alt`, so a low event can be judged against the local
horizon.

`offset_km` is roughly how far the site sits from the transit centerline, and
`path_halfwidth_km` is half the width of the ground strip that sees a transit.
Both are flat-ground estimates from the station's range and altitude, good for
"how close is it" rather than for picking a spot. Use a dedicated transit
site such as nightastro.fr to plan the actual shoot.
