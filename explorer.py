from sympy import symbols, S, diff, limit, oo, solve, re, I, Union, Interval, nan, EmptySet
from sympy.calculus.util import continuous_domain
from sympy.plotting import plot
from sympy.calculus.singularities import singularities
from sympy.parsing.sympy_parser import parse_expr
import sympy
from functools import reduce

x = symbols('x')

class explorer:
    def __init__(self, function):  
        self.function = explorer.get_sym(function)
        self.derivative = self.function.diff(x)
        self.derived_second = self.derivative.diff(x)
        self.domain = [continuous_domain(self.function, x, S.Reals)]
        self.y_point_of_intersection = (0, self.function.subs(x, 0))
        self.x_point_of_intersections = []        
        self.extremums = []
        self.inflection_points = []
        self.vertical_asymptote = []
        self.horizontal_asymptotes = []
        self.place_ud_ = []
        self.place_pn_ = []

        self.x_point_of_intersection()
        self.get_extremum()
        self.get_inflection_point()
        self.find_vertical_asymptote()
        self.find_horizontal_asymptote()
        # למשתנה שבשורה הבאה נוספו ערכי נקודת פיתול מכיוון שיכול להיות שהממוצע יצא אחד מהם
        self.list_of_changes = sorted([-oo, oo] + [xy[0] for m,xy in self.extremums] + self.vertical_asymptote + [xy[0] for un, xy, nu in self.inflection_points] + [xy[0] for xy in self.x_point_of_intersections])
        self.place_pn()
        self.place_ud()
    
    @staticmethod
    def get_sym(str_func):
        new_str_func = ''
        for num_note in range(len(str_func)-1):
            if str_func[num_note] == ')' and str_func[num_note+1]  == '(':
                new_str_func += ')*'
            elif str_func[num_note].isnumeric() and str_func[num_note+1] == 'x':
                new_str_func += str_func[num_note] + '*'
            elif str_func[num_note].isnumeric() and str_func[num_note+1] == '(' or str_func[num_note] == 'x' and str_func[num_note+1] == '(':
                new_str_func += str_func[num_note] + '*'
            else:
                new_str_func += str_func[num_note]
        new_str_func += str_func[-1]   
        return parse_expr(new_str_func)

    def remove_simulated_number(self, lst):
        return list(filter(lambda x: re(x) == x, lst)) 
        
    def x_point_of_intersection(self):
        for solution in self.remove_simulated_number(solve(self.function, x)):
            self.x_point_of_intersections.append((solution, 0))

    def get_extremum(self):
        for suspect in self.remove_simulated_number(solve(self.derivative, x)):
            if self.derivative.diff(x).subs(x, suspect) < 0:
                self.extremums.append(('max', (suspect, self.function.subs(x, suspect))))
            elif self.derivative.diff(x).subs(x, suspect) > 0:
                self.extremums.append(('min', (suspect, self.function.subs(x, suspect))))

    def marking_concavities(self, concavite):
        if concavite == 'u':
            return u'\u22c3'
        return u'\u22c2'

    def get_inflection_point(self):
        list_of_concavites = []
        inflection_points = self.remove_simulated_number(solve(self.derived_second, x)) 
        x_points = sorted([-oo, oo] + inflection_points)
        for i in range(len(x_points)-1):
            average_ = self.between(x_points[i], x_points[i+1])
            if self.derived_second.subs(x, average_) > 0:
                list_of_concavites.append('u')
            elif self.derived_second.subs(x, average_) < 0:
                list_of_concavites.append('n')
        for number, inflection_point in enumerate(inflection_points):
            if list_of_concavites[number] != list_of_concavites[number+1]:
                self.inflection_points.append((self.marking_concavities(list_of_concavites[number]), (inflection_point, self.function.subs(x, inflection_point)), self.marking_concavities(list_of_concavites[number+1])))

    def find_vertical_asymptote(self):
        for non_definition in singularities(self.function, x).args:
            if limit(self.function, x, non_definition) in  [oo, -oo]:
                self.vertical_asymptote.append(non_definition)
            else:
                self.domain.append(u'X\u2260'+str(non_definition))
    
    def find_horizontal_asymptote(self):
        if limit(self.function, x, -oo) not in [-oo, oo]:
            self.horizontal_asymptotes.append('-X -> ' + str(limit(self.function, x, -oo)))
        if limit(self.function, x, oo)  not in [-oo, oo]:
            self.horizontal_asymptotes.append('X -> ' + str(limit(self.function, x, oo)))  

    def between(self, u, v):
        center = (u+v)/2
        if center == oo:
            return min(u, v) + 1
        elif center == -oo:
            return max(u, v) - 1
        return center

    def unify(self, intervals):
        if not intervals:
            return EmptySet
        return reduce(Union, intervals)

    def place_pn(self):  
        list_of_positivity = []
        list_of_negativity = []
        for i in range(len(self.list_of_changes)-1):          
            center = self.between(self.list_of_changes[i], self.list_of_changes[i+1])          
            theorem = self.function.subs(x, center)
            if center == nan:
                if theorem > 0:
                    return ('Permanent positive', EmptySet)
                elif theorem < 0:
                    return (EmptySet, 'Permanent negative')
                elif theorem == 0:
                    return ('Permanent ziro', 'Permanent ziro')
            if theorem > 0:
                list_of_positivity.append(Interval(self.list_of_changes[i], self.list_of_changes[i+1]))
            elif theorem < 0:
                list_of_negativity.append(Interval(self.list_of_changes[i], self.list_of_changes[i+1]))
        self.place_pn_ += [self.unify(list_of_positivity), self.unify(list_of_negativity)] 

    def place_ud(self):  
        list_of_increases = []
        list_of_decreases = []
        for i in range(len(self.list_of_changes)-1):      
            small, large = self.list_of_changes[i], self.list_of_changes[i+1]    
            center = self.between(small, large)  
            theorem = self.derivative.subs(x, center)
            if theorem > 0:
                list_of_increases.append(Interval(small, large))
            elif theorem < 0:
                list_of_decreases.append(Interval(small, large))
        self.place_ud_ += [self.unify(list_of_increases), self.unify(list_of_decreases)]

    def complementary(self):
        complements = [self.function,
            self.derivative,
            self.domain,
            self.y_point_of_intersection,
            self.x_point_of_intersections,
            self.extremums,
            self.vertical_asymptote,            
            self.horizontal_asymptotes,
            self.inflection_points,
            self.place_pn_[0],
            self.place_pn_[1],
            self.place_ud_[0],
            self.place_ud_[1]]
        for num, one in enumerate(complements):
            if one == []:
                complements[num] = None
        return complements

    def __repr__(self):
        return """
            f(x) :{}
            f'(x): {}
            Domain: {}
            Intersection with Y axis : {}
            Intersection with X axis : {}
            Extremums : {}
            Vertical asymptote : {}
            Horizontal asymptote : {}
            inflection_points : {}
            Positive intrvals : {}
            Negative intrvals : {}
            Ascending intervals : {}
            Descending intervals : {}
            """.format(*self.complementary()) 

 
function = '(x**2)/(3(x-2)(x-4))'

print(explorer(function))
plot(explorer.get_sym(function))
