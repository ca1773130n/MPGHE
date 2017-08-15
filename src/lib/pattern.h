#pragma once
#include <string>

namespace utils {
	template <typename T>
		class ISingleton {
			public:
				static T *getInstance(void) {
					if (!mInstance) {
						mInstance = new T();
					}
					return mInstance;
				}

				static void destroyInstance(void) {
					if (mInstance)
						delete mInstance;
				}

			protected:
				ISingleton() {}
				virtual ~ISingleton() {}

			private:
				static T *mInstance;
		};

	template <typename T>
		T * ISingleton<T>::mInstance = 0;
}
