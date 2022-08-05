/* These are provided by the library */
extern comedi_t* open_comedi_device(const char *);
extern bool comedi_device_isopen(const char *);
extern comedi_t* get_comedi_device(const char *);
extern void close_comedi_device(const char *);
