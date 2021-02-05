from sympy import symbols, S, diff, limit, oo, solve, re, I, Union, Interval, nan, EmptySet, sqrt, Reals, zoo
from sympy.calculus.util import continuous_domain
from sympy.plotting import plot
from sympy.calculus.singularities import singularities
from sympy.parsing.sympy_parser import parse_expr
from functools import reduce

x = symbols('x')

class Explorer:
    def __init__(self, function):  
        self.function = Explorer.get_sym(function)
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

    u = u'\u22c2'
    n = u'\u22c3'
    def marking_concavities(self, concavite):
        if concavite == 'u':
            return Explorer.u
        return Explorer.n

    def between(self, num0, num1):
        center = (num0 + num1) / 2
        if center == oo:
            return min(num0, num1) + 1
        elif center == -oo:
            return max(num0, num1) - 1
        elif center == nan: # num0 == -oo and num1 == oo
            return 0
        return center

    def get_inflection_point(self):
        list_of_concavites = []
        inflection_points = self.remove_simulated_number(solve(self.derived_second, x)) 
        if not inflection_points:
            self.inflection_points = []
            return
        x_points = sorted([-oo, oo] + inflection_points )
        for i in range(len(x_points)-1):
            average_ = self.between(x_points[i], x_points[i+1])
            derived2 = self.derived_second.subs(x, average_)
            if derived2 > 0:
                list_of_concavites.append('u')
            elif derived2 < 0:
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
                    return ('Permanent zero', 'Permanent zero')
            if re(theorem) != theorem:
                continue
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
            if re(theorem) != theorem:
                continue
            if theorem > 0:
                list_of_increases.append(Interval(small, large))
            elif theorem < 0:
                list_of_decreases.append(Interval(small, large))
        self.place_ud_ += [self.unify(list_of_increases), self.unify(list_of_decreases)]

    def show_exploration (self):
        inflections = []
        for single_inflection in  self.inflection_points:
            inflections.append((single_inflection[0],[float(a.n()) for a in single_inflection[1]], single_inflection[2]))

        parameters = [self.function,
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
        for num, parameter in enumerate(parameters):
            if parameter == []:
                parameters[num] = None
            elif parameter == Interval(-oo, oo):
                parameters[num] = Reals
        return parameters

    def __repr__(self):
        return """
            f(x) : {}
            f'(x) : {}
            Domain : {}
            Intersection with Y axis : {}
            Intersections with X axis : {}
            Extremums : {}
            Vertical asymptote : {}
            Horizontal asymptote : {}
            inflection_points : {}
            Positive intrvals : {}
            Negative intrvals : {}
            Ascending intervals : {}
            Discending intervals : {}
            """.format(*self.show_exploration()) 

 
# function = 'sqrt(x)'
# a = Explorer(function)
# print(a.show_exploration())
# print(a)
# plot(Explorer.get_sym(function))
