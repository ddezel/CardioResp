/* define useful functions */

#define pi acos(-1.0)

#define max(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a > _b ? _a : _b; })                                        

#define min(a,b)              \
({ __typeof__ (a) _a = (a);   \
  __typeof__ (b) _b = (b);    \
  _a < _b ? _a : _b; })                   