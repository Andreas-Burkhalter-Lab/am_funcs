/** \brief A class that allows sorting an index vector based on a lookup value
 * \tparam CompData A type containing the comparison data
 */
template <typename CompData>
class IndirectLess {
public:
  /** Constructor
   * \param data A vector containing the lookup value for each entry
   */
  IndirectLess(const CompData &data) : mrData(data) {;}
  /** The comparison operator
   * \param i1 The index of the first entry
   * \param i2 The index of the second entry
   * \return True if data[i1] < data[i2]
   */
  bool operator()(int i1,int i2) { return (mrData[i1] < mrData[i2]); }
private:
  const CompData &mrData;
};

/** \brief A class that allows sorting an index vector based on a lookup value, using pointers
 * \tparam CompData A type containing the comparison data
 */
template <typename CompData>
class IndirectLessPtr {
public:
  /** Constructor
   * \param data A vector containing the lookup value for each entry
   */
  IndirectLessPtr(const CompData *data) : mpData(data) {;}
  /** The comparison operator
   * \param i1 The index of the first entry
   * \param i2 The index of the second entry
   * \return True if data[i1] < data[i2]
   */
  bool operator()(int i1,int i2) { return (mpData[i1] < mpData[i2]); }
private:
  const CompData *mpData;
};

/** \brief A class that allows sorting an index vector based on a lookup value
 * \tparam CompData A type containing the comparison data
 */
template <typename CompData>
class IndirectGreater {
public:
  /** Constructor
   * \param data A vector containing the lookup value for each entry
   */
  IndirectGreater(const CompData &data) : mrData(data) {;}
  /** The comparison operator
   * \param i1 The index of the first entry
   * \param i2 The index of the second entry
   * \return True if data[i1] > data[i2]
   */
  bool operator()(int i1,int i2) { return (mrData[i1] > mrData[i2]); }
private:
  const CompData &mrData;
};

/** \brief A class that allows sorting an index vector based on a lookup value, using pointers
 * \tparam CompData A type containing the comparison data
 */
template <typename CompData>
class IndirectGreaterPtr {
public:
  /** Constructor
   * \param data A vector containing the lookup value for each entry
   */
  IndirectGreaterPtr(const CompData *data) : mpData(data) {;}
  /** The comparison operator
   * \param i1 The index of the first entry
   * \param i2 The index of the second entry
   * \return True if data[i1] > data[i2]
   */
  bool operator()(int i1,int i2) { return (mpData[i1] > mpData[i2]); }
private:
  const CompData *mpData;
};
