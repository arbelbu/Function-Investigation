# Function Investigator

This is a python script that investigates mathematical functions.

* Ascending and descending intervals
* Intersections with Y and X axis
* Positive and negative intervals
* Domain
* Extermums
* Asimptotes (vertical and horizontal)
* Derivative

## Example :
```
     input: (x**2)/(3(x-2)(x-4))
     
    output: f(x) : x**2/((x - 4)*(3*x - 6))
            f'(x) : -3*x**2/((x - 4)*(3*x - 6)**2) - x**2/((x - 4)**2*(3*x - 6)) + 2*x/((x - 4)*(3*x - 6))
            Domain : [Union(Interval.open(-oo, 2), Interval.open(2, 4), Interval.open(4, oo))]
            Intersection with Y axis : (0, 0)
            Intersections with X axis : [(0, 0)]
            Extremums : [('min', (0, 0)), ('max', (8/3, -8/3))]
            Vertical asymptote : [2, 4]
            Horizontal asymptote : ['-X -> 1/3', 'X -> 1/3']
            inflection_points : [('⋂', (-4*2**(1/3)/3 - 2*2**(2/3)/3 + 4/3, (-4*2**(1/3)/3 - 2*2**(2/3)/3 + 4/3)**2/((-8/3 - 4*2**(1/3)/3 - 2*2**(2/3)/3)*(-4*2**(1/3) - 2*2**(2/3) - 2))), '⋃')]
            Positive intrvals : Union(Interval(-oo, 2), Interval(4, oo))
            Negative intrvals : Interval(2, 4)
            Ascending intervals : Interval(0, 8/3)
            Descending intervals : Union(Interval(-oo, 0), Interval(8/3, oo))
```
